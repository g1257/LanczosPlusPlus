
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/

#ifndef FEBASED_SC_H
#define FEBASED_SC_H

#include "BasisFeAsBasedSc.h"
#include "SparseRow.h"
#include "ParametersModelFeAs.h"
#include "ModelBase.h"
#include "Geometry/GeometryDca.h"
#include "Parallelizer.h"

namespace LanczosPlusPlus {

template<typename ComplexOrRealType,typename GeometryType,typename InputType>
class FeBasedSc : public ModelBase<ComplexOrRealType,GeometryType,InputType> {

	typedef FeBasedSc<ComplexOrRealType,GeometryType,InputType> ThisType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef ModelBase<ComplexOrRealType,GeometryType,InputType> BaseType;
	typedef PsimagLite::GeometryDca<RealType,GeometryType> GeometryDcaType;

	enum {SPIN_UP = ProgramGlobals::SPIN_UP, SPIN_DOWN = ProgramGlobals::SPIN_DOWN};

	typedef ParametersModelFeAs<RealType> ParametersModelType;
	typedef BasisFeAsBasedSc<GeometryType> BasisType;
	typedef typename BasisType::BaseType BasisBaseType;
	typedef typename BasisType::WordType WordType;
	typedef typename BaseType::SparseMatrixType SparseMatrixType;
	typedef typename BaseType::VectorType VectorType;

	class MatrixVectorHelper {

		typedef PsimagLite::Concurrency ConcurrencyType;

	public:

		MatrixVectorHelper(SizeType nthreads,
		                   VectorType &x,
		                   const VectorType& y,
		                   const BasisBaseType& basis,
		                   const ThisType& myself)
		    : nthreads_(nthreads),x_(x),y_(y),basis_(basis),myself_(myself)
		{}

		void thread_function_(SizeType threadNum,
		                      SizeType blockSize,
		                      SizeType total,
		                      ConcurrencyType::MutexType*)
		{
			SizeType nsite = myself_.geometry_.numberOfSites();
			for (SizeType p=0;p<blockSize;p++) {
				SizeType ispace = threadNum*blockSize + p;
				if (ispace>=total) break;
				SparseRowType sparseRow;

				WordType ket1 = basis_.operator()(ispace,SPIN_UP);
				WordType ket2 = basis_.operator()(ispace,SPIN_DOWN);

				x_[ispace] += myself_.findS(nsite,ket1,ket2,ispace,basis_)*y_[ispace];

				for (SizeType i=0;i<nsite;i++) {
					for (SizeType orb=0;orb<myself_.mp_.orbitals;orb++) {
						myself_.setHoppingTerm(sparseRow,ket1,ket2,
						                       i,orb,basis_);

						if (myself_.mp_.feAsMode == 0) {
							myself_.setU2OffDiagonalTerm(sparseRow,ket1,ket2,
							                             i,orb,basis_);

							for (SizeType orb2=orb+1;orb2<myself_.mp_.orbitals;orb2++) {
								myself_.setU3Term(sparseRow,ket1,ket2,
								                  i,orb,orb2,basis_);
							}

							myself_.setJTermOffDiagonal(sparseRow,ket1,ket2,
							                            i,orb,basis_);
						} else if (myself_.mp_.feAsMode == 1 || myself_.mp_.feAsMode == 2) {
							myself_.setOffDiagonalDecay(sparseRow,ket1,ket2,
							                            i,orb,basis_);
						} else if (myself_.mp_.feAsMode == 3) {
							myself_.setOffDiagonalJimpurity(sparseRow,ket1,ket2,i,orb,basis_);
						} else if (myself_.mp_.feAsMode == 4) {
							myself_.setOffDiagonalKspace(sparseRow,ket1,ket2,i,orb,basis_);
						}
					}
				}

				x_[ispace] += sparseRow.finalize(y_);
			}
		}

	private:

		SizeType nthreads_;
		VectorType &x_;
		const VectorType& y_;
		const BasisBaseType& basis_;
		const ThisType& myself_;
	}; // class MatrixVectorHelper

public:

	typedef PsimagLite::SparseRow<SparseMatrixType> SparseRowType;
	enum {TERM_HOPPINGS = 0,TERM_J_PM = 1, TERM_J_ZZ = 2};
	static int const FERMION_SIGN = BasisType::FERMION_SIGN;

	FeBasedSc(SizeType nup,SizeType ndown,InputType& io,const GeometryType& geometry)
	    : mp_(io),
	      geometry_(geometry),
	      basis_(geometry,nup,ndown,mp_.orbitals),
	      geometryDca_(geometry,mp_.orbitals)
	{}

	~FeBasedSc()
	{
		BaseType::deleteGarbage(garbage_);
	}

	SizeType size() const { return basis_.size(); }

	SizeType orbitals(SizeType) const
	{
		return mp_.orbitals;
	}

	void setupHamiltonian(SparseMatrixType& matrix) const
	{
		setupHamiltonian(matrix,basis_);
	}

	bool hasNewParts(std::pair<SizeType,SizeType>& newParts,
	                 SizeType what,
	                 SizeType spin,
	                 const std::pair<SizeType,SizeType>& orbs) const
	{
		return basis_.hasNewParts(newParts,what,spin,orbs);
	}

	const GeometryType& geometry() const { return geometry_; }

	void setupHamiltonian(SparseMatrixType& matrix,
	                      const BasisBaseType& basis) const
	{
		// Calculate diagonal elements AND count non-zero matrix elements
		SizeType hilbert=basis.size();
		typename PsimagLite::Vector<RealType>::Type diag(hilbert);
		calcDiagonalElements(diag,basis);

		SizeType nsite = geometry_.numberOfSites();

		// Setup CRS matrix
		matrix.resize(hilbert,hilbert);

		// Calculate off-diagonal elements AND store matrix
		SizeType nCounter=0;
		for (SizeType ispace=0;ispace<hilbert;ispace++) {
			SparseRowType sparseRow;
			matrix.setRow(ispace,nCounter);
			WordType ket1 = basis(ispace,SPIN_UP);
			WordType ket2 = basis(ispace,SPIN_DOWN);
			// Save diagonal
			sparseRow.add(ispace,diag[ispace]);
			for (SizeType i=0;i<nsite;i++) {
				for (SizeType orb=0;orb<mp_.orbitals;orb++) {
					setHoppingTerm(sparseRow,ket1,ket2,i,orb,basis);

					if (mp_.feAsMode == 0) {
						setU2OffDiagonalTerm(sparseRow,ket1,ket2,
						                     i,orb,basis);
						for (SizeType orb2=0;orb2<mp_.orbitals;orb2++) {
							if (orb==orb2) continue;

							setU3Term(sparseRow,ket1,ket2,
							          i,orb,orb2,basis);
						}

						setJTermOffDiagonal(sparseRow,ket1,ket2,
						                    i,orb,basis);
					} else if (mp_.feAsMode == 1 || mp_.feAsMode == 2) {
						setOffDiagonalDecay(sparseRow,ket1,ket2,
						                    i,orb,basis);
					} else if (mp_.feAsMode == 3) {
						setOffDiagonalJimpurity(sparseRow,ket1,ket2,i,orb,basis);
					} else if (mp_.feAsMode == 4) {
						setOffDiagonalKspace(sparseRow,ket1,ket2,i,orb,basis);
					}
				}
			}

			nCounter += sparseRow.finalize(matrix);
		}

		matrix.setRow(hilbert,nCounter);
	}

	void matrixVectorProduct(VectorType &x,const VectorType& y) const
	{
		matrixVectorProduct(x,y,basis_);
	}

	void matrixVectorProduct(VectorType &x,
	                         const VectorType& y,
	                         const BasisBaseType& basis) const
	{
		SizeType hilbert=basis.size();

		// Calculate off-diagonal elements AND store matrix
		typedef MatrixVectorHelper HelperType;
		typedef PsimagLite::Parallelizer<HelperType> ParallelizerType;
		ParallelizerType threadObject(PsimagLite::Concurrency::npthreads,
		                              PsimagLite::MPI::COMM_WORLD);
		HelperType helper(PsimagLite::Concurrency::npthreads,x,y,basis,*this);

		std::cout<<"Using "<<threadObject.name();
		std::cout<<" with "<<threadObject.threads()<<" threads.\n";
		threadObject.loopCreate(hilbert,helper);
	}

	const BasisType& basis() const { return basis_; }

	PsimagLite::String name() const { return __FILE__; }

	BasisType* createBasis(SizeType nup, SizeType ndown) const
	{
		BasisType* ptr = new BasisType(geometry_,nup,ndown);
		garbage_.push_back(ptr);
		return ptr;
	}

	void print(std::ostream& os) const { os<<mp_; }

private:

	void setOffDiagonalDecay(SparseRowType& sparseRow,
	                         const WordType& ket1,
	                         const WordType& ket2,
	                         SizeType i,
	                         SizeType,
	                         const BasisBaseType &basis) const
	{
		for (SizeType spin1 = 0; spin1 < 2; ++spin1) {
			for (SizeType spin2 = 0; spin2 < 2; ++spin2) {

				if (spin1 != spin2) continue;

				if (!basis.isThereAnElectronAt(ket1,ket2,i,spin1,1)) continue;
				if (!basis.isThereAnElectronAt(ket1,ket2,i,spin2,1)) continue;
				WordType mask = BasisType::bitmask(i*mp_.orbitals+1);
				WordType bra1 = (spin1 == SPIN_UP) ? (ket1 ^ mask): ket1;
				WordType bra2 = (spin1 == SPIN_UP) ? ket2 : (ket2 ^ mask);

				if (!basis.isThereAnElectronAt(bra1,bra2,i,spin2,1)) continue;
				WordType bra3 = (spin2 == SPIN_UP) ? (bra1 ^ mask): bra1;
				WordType bra4 = (spin2 == SPIN_UP) ? bra2 : (bra2 ^ mask);

				if (basis.isThereAnElectronAt(bra3,bra4,i,spin2,0)) continue;
				mask = BasisType::bitmask(i*mp_.orbitals+0);
				WordType bra5 = (spin2 == SPIN_UP) ? (bra3 ^ mask): bra3;
				WordType bra6 = (spin2 == SPIN_UP) ? bra4 : (bra4 ^ mask);

				if (basis.isThereAnElectronAt(bra5,bra6,i,spin1,2)) continue;
				mask = BasisType::bitmask(i*mp_.orbitals+2);
				WordType bra7 = (spin1 == SPIN_UP) ? (bra5 ^ mask): bra5;
				WordType bra8 = (spin1 == SPIN_UP) ? bra6 : (bra6 ^ mask);

				SizeType temp = basis.perfectIndex(bra7,bra8);
				sparseRow.add(temp,mp_.coulombV);
			}
		}
	}

	RealType findSdecay(SizeType,
	                    WordType ket1,
	                    WordType ket2,
	                    SizeType i,
	                    SizeType orb,
	                    const BasisBaseType& basis) const
	{
		RealType s = mp_.hubbardU[orb+orb*mp_.orbitals] *
		        basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb) *
		        basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb);

		for (SizeType orb2=orb+1;orb2<mp_.orbitals;orb2++) {
			// Hubbard term U1
			s += mp_.hubbardU[orb+orb2*mp_.orbitals] * nix(ket1,ket2,i,orb,basis) *
			        nix(ket1,ket2,i,orb2,basis);
		}

		return s;
	}

	ComplexOrRealType hoppings(SizeType i,SizeType orb1,SizeType j,SizeType orb2) const
	{
		return -geometry_(i,orb1,j,orb2,TERM_HOPPINGS);
	}

	void setHoppingTerm(SparseRowType &sparseRow,
	                    const WordType& ket1,
	                    const WordType& ket2,
	                    SizeType i,
	                    SizeType orb,
	                    const BasisBaseType &basis) const
	{
		SizeType ii = i*mp_.orbitals+orb;
		WordType s1i=(ket1 & BasisType::bitmask(ii));
		if (s1i>0) s1i=1;
		WordType s2i=(ket2 & BasisType::bitmask(ii));
		if (s2i>0) s2i=1;

		SizeType nsite = geometry_.numberOfSites();

		// Hopping term
		for (SizeType j=0;j<nsite;j++) {
			if (j<i) continue;
			for (SizeType orb2=0;orb2<mp_.orbitals;orb2++) {
				SizeType jj = j*mp_.orbitals+orb2;
				ComplexOrRealType h = hoppings(i,orb,j,orb2);
				if (PsimagLite::real(h) == 0 && PsimagLite::imag(h) == 0) continue;
				WordType s1j= (ket1 & BasisType::bitmask(jj));
				if (s1j>0) s1j=1;
				WordType s2j= (ket2 & BasisType::bitmask(jj));
				if (s2j>0) s2j=1;

				if (s1i+s1j==1) {
					WordType bra1= ket1 ^(BasisType::bitmask(ii)|BasisType::bitmask(jj));
					SizeType temp = basis.perfectIndex(bra1,ket2);
					RealType extraSign = (s1i==1) ? FERMION_SIGN : 1;
					RealType tmp2 = basis_.doSign(ket1,ket2,i,orb,j,orb2,SPIN_UP);
					ComplexOrRealType cTemp = h*extraSign*tmp2;
					sparseRow.add(temp,cTemp);

				}
				if (s2i+s2j==1) {
					WordType bra2= ket2 ^(BasisType::bitmask(ii)|BasisType::bitmask(jj));
					SizeType temp = basis.perfectIndex(ket1,bra2);
					RealType extraSign = (s2i==1) ? FERMION_SIGN : 1;
					RealType tmp2 = basis_.doSign(ket1,ket2,i,orb,j,orb2,SPIN_DOWN);
					ComplexOrRealType cTemp = h*extraSign*tmp2;
					sparseRow.add(temp,cTemp);
				}
			}
		}
	}

	void setU2OffDiagonalTerm(
	        SparseRowType &sparseRow,
	        const WordType& ket1,
	        const WordType& ket2,
	        SizeType i,
	        SizeType orb1,
	        const BasisBaseType &basis) const
	{
		RealType val = FERMION_SIGN * mp_.hubbardU[2]*0.5;

		for (SizeType orb2=orb1+1;orb2<mp_.orbitals;orb2++) {
			setSplusSminus(sparseRow,ket1,ket2,i,orb1,i,orb2,val,basis);
			setSplusSminus(sparseRow,ket1,ket2,i,orb2,i,orb1,val,basis);
		}

	}

	// N.B.: orb1!=orb2 here
	void setSplusSminus(
	        SparseRowType &sparseRow,
	        const WordType& ket1,
	        const WordType& ket2,
	        SizeType i,
	        SizeType orb1,
	        SizeType j,
	        SizeType orb2,
	        ComplexOrRealType value,
	        const BasisBaseType &basis) const
	{
		if (splusSminusNonZero(ket1,ket2,i,orb1,j,orb2,basis)==0) return;

		SizeType ii = i*mp_.orbitals + orb1;
		SizeType jj = j*mp_.orbitals + orb2;
		assert(ii!=jj);
		WordType bra1 = ket1 ^ (BasisType::bitmask(ii)|BasisType::bitmask(jj));
		WordType bra2 = ket2 ^ (BasisType::bitmask(ii)|BasisType::bitmask(jj));
		SizeType temp = basis.perfectIndex(bra1,bra2);
		sparseRow.add(temp,value);
	}

	// N.B.: orb1!=orb2 here
	void setU3Term(
	        SparseRowType &sparseRow,
	        const WordType& ket1,
	        const WordType& ket2,
	        SizeType i,
	        SizeType orb1,
	        SizeType orb2,
	        const BasisBaseType &basis) const
	{
		assert(orb1!=orb2);
		if (u3TermNonZero(ket1,ket2,i,orb1,orb2,basis)==0) return;

		SizeType ii = i*mp_.orbitals + orb1;
		SizeType jj = i*mp_.orbitals + orb2;
		WordType bra1 = ket1 ^ (BasisType::bitmask(ii)|BasisType::bitmask(jj));
		WordType bra2 = ket2 ^ (BasisType::bitmask(ii)|BasisType::bitmask(jj));
		SizeType temp = basis.perfectIndex(bra1,bra2);
		sparseRow.add(temp,FERMION_SIGN * mp_.hubbardU[3]);
	}

	void setJTermOffDiagonal(
	        SparseRowType &sparseRow,
	        const WordType& ket1,
	        const WordType& ket2,
	        SizeType i,
	        SizeType orb,
	        const BasisBaseType& basis) const
	{
		for (SizeType j=0;j<geometry_.numberOfSites();j++) {
			ComplexOrRealType value = jCoupling(i,j,TERM_J_PM)*0.5;
			if (PsimagLite::real(value) == 0 && PsimagLite::imag(value) == 0) continue;
			value *= 0.5; // RealType counting i,j
			assert(i!=j);
			for (SizeType orb2=0;orb2<mp_.orbitals;orb2++) {
				//if (orb2!=orb) continue; // testing only!!
				RealType sign = jTermSign(ket1,ket2,i,orb,j,orb2,basis);
				setSplusSminus(sparseRow,ket1,ket2,i,orb,j,orb2,value*sign,basis);
				setSplusSminus(sparseRow,ket1,ket2,j,orb2,i,orb,value*sign,basis);
			}
		}
	}

	int jTermSign(
	        const WordType& ket1,
	        const WordType& ket2,
	        SizeType i,
	        SizeType orb1,
	        SizeType j,
	        SizeType orb2,
	        const BasisBaseType &basis) const
	{
		if (i>j) return jTermSign(ket1,ket2,j,orb2,i,orb1,basis);
		int x = basis.doSign(ket1,ket2,i,orb1,j,orb2,SPIN_UP);
		x *= basis.doSign(ket1,ket2,i,orb1,j,orb2,SPIN_DOWN);
		return x;
	}

	void calcDiagonalElements(typename PsimagLite::Vector<RealType>::Type& diag,
	                          const BasisBaseType& basis) const
	{
		SizeType hilbert=basis.size();
		SizeType nsite = geometry_.numberOfSites();

		// Calculate diagonal elements
		for (SizeType ispace=0;ispace<hilbert;ispace++) {
			WordType ket1 = basis(ispace,SPIN_UP);
			WordType ket2 = basis(ispace,SPIN_DOWN);
			diag[ispace]=findS(nsite,ket1,ket2,ispace,basis);
		}
	}

	RealType findS(SizeType nsite,
	               WordType ket1,
	               WordType ket2,
	               SizeType,
	               const BasisBaseType& basis) const
	{
		RealType s = 0;
		for (SizeType i=0;i<nsite;i++) {
			for (SizeType orb=0;orb<mp_.orbitals;orb++) {

				if (mp_.feAsMode == 0) {
					s += findSnoDecay(nsite,ket1,ket2,i,orb,basis);
				} else if (mp_.feAsMode == 1 || mp_.feAsMode == 2){
					s += findSdecay(nsite,ket1,ket2,i,orb,basis);
				} else if (mp_.feAsMode == 3) {
					s += findSImpurity(nsite,ket1,ket2,i,orb,basis);
				} else if (mp_.feAsMode == 4) {
					s += findSkspace(nsite,ket1,ket2,i,orb,basis);
				}

				// Potential term
				s += mp_.potentialV[i+(orb+mp_.orbitals*0)*nsite]*
				        basis.getN(ket1,ket1,i,SPIN_UP,orb) +
				        mp_.potentialV[i+(orb+mp_.orbitals*1)*nsite]*
				        basis.getN(ket2,ket2,i,SPIN_DOWN,orb);

			}
		}
		return s;
	}

	RealType findSnoDecay(SizeType nsite,
	                      WordType ket1,
	                      WordType ket2,
	                      SizeType i,
	                      SizeType orb,
	                      const BasisBaseType& basis) const
	{
		// Hubbard term U0
		ComplexOrRealType s = mp_.hubbardU[0]*basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb)
		        * basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb);

		for (SizeType orb2=orb+1;orb2<mp_.orbitals;orb2++) {
			// Hubbard term U1
			s += mp_.hubbardU[1] * nix(ket1,ket2,i,orb,basis) *
			        nix(ket1,ket2,i,orb2,basis);

			// Diagonal U2 term
			s+= mp_.hubbardU[4]*
			        szTerm(ket1,ket2,i,orb,basis)*
			        szTerm(ket1,ket2,i,orb2,basis);
		}

		// JNN and JNNN diagonal part
		for (SizeType j=0;j<nsite;j++) {
			for (SizeType orb2=0;orb2<mp_.orbitals;orb2++) {
				ComplexOrRealType value = jCoupling(i,j,TERM_J_ZZ);
				if (PsimagLite::real(value) == 0 && PsimagLite::imag(value) == 0) continue;
				s += value*0.5* // RealType counting i,j
				        szTerm(ket1,ket2,i,orb,basis)*
				        szTerm(ket1,ket2,j,orb2,basis);
			}
		}

		assert(fabs(PsimagLite::imag(s))<1e-12);
		return PsimagLite::real(s);
	}

	RealType findSImpurity(SizeType,
	                       WordType ket1,
	                       WordType ket2,
	                       SizeType i,
	                       SizeType orb,
	                       const BasisBaseType& basis) const
	{
		if (i > 0) return 0.0;

		// Hubbard term U0
		RealType s = mp_.hubbardU[0]*basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb)
		        *basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb);

		for (SizeType orb2=0;orb2<mp_.orbitals;orb2++) {
			// Hubbard term U1
			if (orb == orb2) continue;
			for (SizeType spin = 0; spin < 2; ++spin)
				s += 0.5*mp_.hubbardU[1]*basis.isThereAnElectronAt(ket1,ket2,i,spin,orb)
				        *basis.isThereAnElectronAt(ket1,ket2,i,spin,orb2);
		}

		for (SizeType orb2=0;orb2<mp_.orbitals;orb2++) {
			if (orb == orb2) continue;
			// Diagonal U2 term
			s+= mp_.hubbardU[4]*
			        basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb)
			        *basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb2);
		}

		return s;
	}

	RealType findSkspace(SizeType,
	                     WordType ket1,
	                     WordType ket2,
	                     SizeType i,
	                     SizeType orb,
	                     const BasisBaseType& basis) const
	{
		if (i > 0) return 0.0;
		RealType s = 0;
		SizeType ck1 = basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb);
		if (ck1 == 0) return 0.0;

		for (SizeType orb2=0;orb2<mp_.orbitals;orb2++) {
			SizeType ck2 = basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb2);
			if (ck2 == 0) continue;
			s++;
		}

		return s*mp_.hubbardU[0];
	}

	SizeType splusSminusNonZero(
	        const WordType& ket1,
	        const WordType& ket2,
	        SizeType i,
	        SizeType orb1,
	        SizeType j,
	        SizeType orb2,
	        const BasisBaseType& basis) const
	{
		if (basis.isThereAnElectronAt(ket1,ket2,
		                              j,SPIN_UP,orb2)==0) return 0;
		if (basis.isThereAnElectronAt(ket1,ket2,
		                              i,SPIN_UP,orb1)==1) return 0;
		if (basis.isThereAnElectronAt(ket1,ket2,
		                              i,SPIN_DOWN,orb1)==0) return 0;
		if (basis.isThereAnElectronAt(ket1,ket2,
		                              j,SPIN_DOWN,orb2)==1) return 0;
		return 1;
	}

	SizeType u3TermNonZero(
	        const WordType& ket1,
	        const WordType& ket2,
	        SizeType i,
	        SizeType orb1,
	        SizeType orb2,
	        const BasisBaseType &basis) const
	{
		if (basis.isThereAnElectronAt(ket1,ket2,
		                              i,SPIN_UP,orb2)==0) return 0;
		if (basis.isThereAnElectronAt(ket1,ket2,
		                              i,SPIN_UP,orb1)==1) return 0;
		if (basis.isThereAnElectronAt(ket1,ket2,
		                              i,SPIN_DOWN,orb1)==1) return 0;
		if (basis.isThereAnElectronAt(ket1,ket2,
		                              i,SPIN_DOWN,orb2)==0) return 0;
		return 1;
	}

	SizeType nix(
	        const WordType& ket1,
	        const WordType& ket2,
	        SizeType i,
	        SizeType orb,
	        const BasisBaseType& basis) const
	{
		SizeType sum = 0;
		for (SizeType spin=0;spin<2;spin++)
			sum += basis.isThereAnElectronAt(ket1,ket2,i,spin,orb);
		return sum;
	}

	RealType szTerm(
	        const WordType& ket1,
	        const WordType& ket2,
	        SizeType i,
	        SizeType orb,
	        const BasisBaseType& basis) const
	{
		RealType sz = basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb);
		sz -= basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb);
		return 0.5*sz;
	}

	ComplexOrRealType jCoupling(SizeType i,SizeType j, SizeType term) const
	{
		if (geometry_.terms()==1) return 0.0;
		return geometry_(i,0,j,0,term);
	}

	void setOffDiagonalJimpurity(SparseRowType& sparseRow,
	                             const WordType& ket1,
	                             const WordType& ket2,
	                             SizeType i,
	                             SizeType orb1,
	                             const BasisBaseType &basis) const
	{
		if (i > 0) return;

		SizeType orbitals = mp_.orbitals;

		for (SizeType type = 0; type < 2; ++type) {
			for (SizeType orb2 = 0; orb2 < orbitals; ++orb2) {

				if (orb1 == orb2) continue;

				SizeType orb3 = (type == 0) ? orb2 : orb1;
				SizeType orb4 = (type == 0) ? orb1 : orb2;

				if (!basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb4)) continue;
				if (basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb3)) continue;
				if (!basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb2)) continue;
				if (basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb1)) continue;

				WordType mask4 = BasisType::bitmask(i*mp_.orbitals+orb4);
				WordType mask3 = BasisType::bitmask(i*mp_.orbitals+orb3);
				WordType bra2 =(ket2 ^ mask4) ^ mask3;

				WordType mask2 = BasisType::bitmask(i*mp_.orbitals+orb2);
				WordType mask1 = BasisType::bitmask(i*mp_.orbitals+orb1);
				WordType bra1 = (ket1 ^ mask2) ^ mask1;

				RealType x = basis.doSign(ket1,ket2,i,orb1,i,orb2,SPIN_UP);
				x *= basis.doSign(ket1,ket2,i,orb3,i,orb4,SPIN_DOWN);

				SizeType temp = basis.perfectIndex(bra1,bra2);
				sparseRow.add(temp,x * mp_.hubbardU[3]);
			}
		}
	}

	void setOffDiagonalKspace(SparseRowType& sparseRow,
	                          const WordType& ket1,
	                          const WordType& ket2,
	                          SizeType i,
	                          SizeType orb1,
	                          const BasisBaseType &basis) const
	{
		if (i > 0) return;

		SizeType orbitals = mp_.orbitals;

		if (basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb1)) return;

		for (SizeType orb2 = 0; orb2 < orbitals; ++orb2) {
			if (orb1 == orb2) continue;

			if (!basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb2)) continue;

			for (SizeType orb3 = 0; orb3 < orbitals; ++orb3) {

				if (basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb3)) continue;

				SizeType orb4 = getMomentum(orb1, orb2, orb3);
				assert(orb3 != orb4);

				if (!basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb4)) continue;

				WordType mask4 = BasisType::bitmask(i*mp_.orbitals+orb4);
				WordType mask3 = BasisType::bitmask(i*mp_.orbitals+orb3);
				WordType bra2 =(ket2 ^ mask4) ^ mask3;

				WordType mask2 = BasisType::bitmask(i*mp_.orbitals+orb2);
				WordType mask1 = BasisType::bitmask(i*mp_.orbitals+orb1);
				WordType bra1 = (ket1 ^ mask2) ^ mask1;

				RealType x = basis.doSign(ket1,ket2,i,orb1,i,orb2,SPIN_UP);
				x *= basis.doSign(ket1,ket2,i,orb3,i,orb4,SPIN_DOWN);

				SizeType temp = basis.perfectIndex(bra1,bra2);
				sparseRow.add(temp,x * mp_.hubbardU[0]);
			}
		}
	}

	SizeType getMomentum(SizeType orb1, SizeType orb2, SizeType orb3) const
	{
		assert(orb1 < mp_.orbitals);
		assert(orb2 < mp_.orbitals);
		assert(orb3 < mp_.orbitals);

		SizeType tmp = geometryDca_.kSum(orb3,orb1);
		SizeType orb4 = geometryDca_.kSustract(tmp,orb2);

		assert(orb4 < mp_.orbitals);
		return orb4;
	}

	const ParametersModelType mp_;
	const GeometryType& geometry_;
	BasisType basis_;
	GeometryDcaType geometryDca_;
	mutable typename PsimagLite::Vector<BasisType*>::Type garbage_;
}; // class FeBasedSc

} // namespace LanczosPlusPlus
#endif

