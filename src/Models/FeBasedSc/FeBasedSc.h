
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

#include "CrsMatrix.h"
#include "BasisFeAsBasedSc.h"
#include "SparseRow.h"
#include "ParametersModelFeAs.h"

namespace LanczosPlusPlus {
	
	template<typename RealType_,typename GeometryType_>
	class FeBasedSc {
		
		typedef PsimagLite::Matrix<RealType_> MatrixType;

	public:

		typedef ParametersModelFeAs<RealType_> ParametersModelType;
		typedef GeometryType_ GeometryType;
		typedef BasisFeAsBasedSc<GeometryType> BasisType;
		typedef typename BasisType::WordType WordType;
		typedef RealType_ RealType;
		typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;
		typedef PsimagLite::SparseRow<SparseMatrixType> SparseRowType;
		typedef std::vector<RealType> VectorType;
		enum {SPIN_UP=BasisType::SPIN_UP,SPIN_DOWN=BasisType::SPIN_DOWN};
		enum {DESTRUCTOR=BasisType::DESTRUCTOR,CONSTRUCTOR=BasisType::CONSTRUCTOR};
		enum {TERM_HOPPINGS=0,TERM_J=1};
		static int const FERMION_SIGN = BasisType::FERMION_SIGN;
		
		
		FeBasedSc(size_t nup,size_t ndown,const ParametersModelType& mp,const GeometryType& geometry)
		: mp_(mp),geometry_(geometry),basis_(geometry,nup,ndown,mp_.orbitals)
		{
		}
		
		size_t size() const { return basis_.size(); }

		size_t orbitals(size_t site) const
		{
			return mp_.orbitals;
		}

		void setupHamiltonian(SparseMatrixType &matrix) const
		{
			setupHamiltonian(matrix,basis_);
		}
		
		bool hasNewParts(std::pair<size_t,size_t>& newParts,
						 size_t what2,
						 size_t type,
						 size_t spin,
						 const std::pair<size_t,size_t>& orbs) const
		{
			return basis_.hasNewParts(newParts,type,spin,orbs);
		}

		template<typename SomeVectorType>
		void getModifiedState(SomeVectorType& modifVector,
							  size_t what2,
							  const SomeVectorType& gsVector,
							  const BasisType& basisNew,
							  size_t type,
							  size_t isite,
							  size_t jsite,
							  size_t spin) const
		{
			std::string s = "FeBasedSc::getModifiedState(...): unimplemented. ";
			s+= "This probably means that you can't compute the Green function";
			s+= " with this model (sorry). It might be added in the future.\n";
			throw std::runtime_error(s.c_str());
		}

		const GeometryType& geometry() const { return geometry_; }

		void setupHamiltonian(SparseMatrixType &matrix,
				      const BasisType &basis) const
		{
			// Calculate diagonal elements AND count non-zero matrix elements
			size_t hilbert=basis.size();
			std::vector<RealType> diag(hilbert);
			calcDiagonalElements(diag,basis);

			size_t nsite = geometry_.numberOfSites();

			// Setup CRS matrix
			matrix.resize(hilbert,hilbert);

			// Calculate off-diagonal elements AND store matrix
			size_t nCounter=0;
			for (size_t ispace=0;ispace<hilbert;ispace++) {
				SparseRowType sparseRow;
				matrix.setRow(ispace,nCounter);
				WordType ket1 = basis(ispace,SPIN_UP);
				WordType ket2 = basis(ispace,SPIN_DOWN);
				// Save diagonal
				sparseRow.add(ispace,diag[ispace]);
				for (size_t i=0;i<nsite;i++) {
					for (size_t orb=0;orb<mp_.orbitals;orb++) {
						setHoppingTerm(sparseRow,ket1,ket2,
								i,orb,basis);
						if (orb==0) {
							setU2OffDiagonalTerm(sparseRow,ket1,ket2,
								i,orb,basis);
						}
						setU3Term(sparseRow,ket1,ket2,
								i,orb,1-orb,basis);
						setJTermOffDiagonal(sparseRow,ket1,ket2,
								i,orb,basis);
					}
				}
				nCounter += sparseRow.finalize(matrix);
			}
			matrix.setRow(hilbert,nCounter);
		}

		void matrixVectorProduct(VectorType &x,const VectorType& y) const
		{
			matrixVectorProduct(x,y,&basis_);
		}

		void matrixVectorProduct(VectorType &x,
					  const VectorType& y,
					  const BasisType* basis) const
		{
			// Calculate diagonal elements AND count non-zero matrix elements
			size_t hilbert=basis->size();
			calcDiagonalElements(x,y,basis);

			size_t nsite = geometry_.numberOfSites();

			// Calculate off-diagonal elements AND store matrix

			for (size_t ispace=0;ispace<hilbert;ispace++) {
				SparseRowType sparseRow;

				WordType ket1 = basis->operator ()(ispace,SPIN_UP);
				WordType ket2 = basis->operator ()(ispace,SPIN_DOWN);

				//x[ispace] += diag[ispace]*y[ispace];
				for (size_t i=0;i<nsite;i++) {
					for (size_t orb=0;orb<mp_.orbitals;orb++) {
						setHoppingTerm(sparseRow,ket1,ket2,
								i,orb,*basis);

						setU2OffDiagonalTerm(sparseRow,ket1,ket2,
								i,orb,*basis);

						for (size_t orb2=orb+1;orb2<mp_.orbitals;orb2++) {
							setU3Term(sparseRow,ket1,ket2,
								i,orb,orb2,*basis);
						}
						setJTermOffDiagonal(sparseRow,ket1,ket2,
								i,orb,*basis);
					}
				}
				x[ispace] += sparseRow.finalize(y);
			}
		}

		const BasisType& basis() const { return basis_; }

		template<typename SomeVectorType>
		void accModifiedState(SomeVectorType& z,
							  size_t what2,
							  const BasisType& newBasis,
							  const SomeVectorType& gsVector,
							  size_t what,
							  size_t site,
							  size_t spin,
							  size_t orb,
							  int isign) const
		{
			for (size_t ispace=0;ispace<basis_.size();ispace++) {
				WordType ket1 = basis_(ispace,SPIN_UP);
				WordType ket2 = basis_(ispace,SPIN_DOWN);
				int temp = newBasis.getBraIndex(ket1,ket2,what,site,spin,orb);
				// 				int temp= getBraIndex(mysign,ket1,ket2,newBasis,what,site,spin);
				if (temp>=0 && size_t(temp)>=z.size()) {
					std::string s = "old basis=" + ttos(basis_.size());
					s += " newbasis=" + ttos(newBasis.size());
					s += "\n";
					s += "what=" + ttos(what) + " spin=" + ttos(spin);
					s += " site=" + ttos(site);
					s += "ket1=" + ttos(ket1) + " and ket2=" + ttos(ket2);
					s += "\n";
					s += "getModifiedState: z.size=" + ttos(z.size());
					s += " but temp=" + ttos(temp) + "\n";
					throw std::runtime_error(s.c_str());
				}
				if (temp<0) continue;
//				int mysign = basis_.doSignGf(ket1,ket2,site,spin,orb);
				int mysign = (ProgramGlobals::isFermionic(what2)) ? basis_.doSignGf(ket1,ket2,site,spin,orb) : 1;
				z[temp] += isign*mysign*gsVector[ispace];
			}
		}

	private:

		RealType hoppings(size_t i,size_t orb1,size_t j,size_t orb2) const
		{
			return -geometry_(i,orb1,j,orb2,TERM_HOPPINGS);
		}

		void setHoppingTerm(
				SparseRowType &sparseRow,
				const WordType& ket1,
				const WordType& ket2,
				size_t i,
				size_t orb,
				const BasisType &basis) const
		{
			size_t ii = i*mp_.orbitals+orb;
			WordType s1i=(ket1 & BasisType::bitmask(ii));
			if (s1i>0) s1i=1;
			WordType s2i=(ket2 & BasisType::bitmask(ii));
			if (s2i>0) s2i=1;

			size_t nsite = geometry_.numberOfSites();

			// Hopping term
			for (size_t j=0;j<nsite;j++) {
				if (j<i) continue;
				for (size_t orb2=0;orb2<mp_.orbitals;orb2++) {
					size_t jj = j*mp_.orbitals+orb2;
					RealType h = hoppings(i,orb,j,orb2);
					if (h==0) continue;
					WordType s1j= (ket1 & BasisType::bitmask(jj));
					if (s1j>0) s1j=1;
					WordType s2j= (ket2 & BasisType::bitmask(jj));
					if (s2j>0) s2j=1;

					if (s1i+s1j==1) {
						WordType bra1= ket1 ^(BasisType::bitmask(ii)|BasisType::bitmask(jj));
						size_t temp = basis.perfectIndex(bra1,ket2);
						int extraSign = (s1i==1) ? FERMION_SIGN : 1;
						RealType cTemp = h*extraSign*basis_.doSign(
							ket1,ket2,i,orb,j,orb2,SPIN_UP);
						sparseRow.add(temp,cTemp);

					}
					if (s2i+s2j==1) {
						WordType bra2= ket2 ^(BasisType::bitmask(ii)|BasisType::bitmask(jj));
						size_t temp = basis.perfectIndex(ket1,bra2);
						int extraSign = (s2i==1) ? FERMION_SIGN : 1;
						RealType cTemp = h*extraSign*basis_.doSign(
							ket1,ket2,i,orb,j,orb2,SPIN_DOWN);
						sparseRow.add(temp,cTemp);
					}
				}
			}
		}
		
		void setU2OffDiagonalTerm(
				SparseRowType &sparseRow,
				const WordType& ket1,
				const WordType& ket2,
				size_t i,
				size_t orb1,
				const BasisType &basis) const
		{
			RealType val = FERMION_SIGN * mp_.hubbardU[2]*0.5;

			for (size_t orb2=orb1+1;orb2<mp_.orbitals;orb2++) {
				setSplusSminus(sparseRow,ket1,ket2,i,orb1,i,orb2,val,basis);
				setSplusSminus(sparseRow,ket1,ket2,i,orb2,i,orb1,val,basis);
			}

		}

		// N.B.: orb1!=orb2 here
		void setSplusSminus(
				SparseRowType &sparseRow,
				const WordType& ket1,
				const WordType& ket2,
				size_t i,
				size_t orb1,
				size_t j,
				size_t orb2,
				RealType value,
				const BasisType &basis) const
		{
			if (splusSminusNonZero(ket1,ket2,i,orb1,j,orb2,basis)==0) return;

			size_t ii = i*mp_.orbitals + orb1;
			size_t jj = j*mp_.orbitals + orb2;
			assert(ii!=jj);
			WordType bra1 = ket1 ^ (BasisType::bitmask(ii)|BasisType::bitmask(jj));
			WordType bra2 = ket2 ^ (BasisType::bitmask(ii)|BasisType::bitmask(jj));
			size_t temp = basis.perfectIndex(bra1,bra2);
			sparseRow.add(temp,value);
		}

		// N.B.: orb1!=orb2 here
		void setU3Term(
				SparseRowType &sparseRow,
				const WordType& ket1,
				const WordType& ket2,
				size_t i,
				size_t orb1,
				size_t orb2,
				const BasisType &basis) const
		{
			assert(orb1!=orb2);
			if (u3TermNonZero(ket1,ket2,i,orb1,orb2,basis)==0) return;

			size_t ii = i*mp_.orbitals + orb1;
			size_t jj = i*mp_.orbitals + orb2;
			WordType bra1 = ket1 ^ (BasisType::bitmask(ii)|BasisType::bitmask(jj));
			WordType bra2 = ket2 ^ (BasisType::bitmask(ii)|BasisType::bitmask(jj));
			size_t temp = basis.perfectIndex(bra1,bra2);
			sparseRow.add(temp,FERMION_SIGN * mp_.hubbardU[3]);
		}

		void setJTermOffDiagonal(
				SparseRowType &sparseRow,
				const WordType& ket1,
				const WordType& ket2,
				size_t i,
				size_t orb,
				const BasisType &basis) const
		{
			for (size_t j=0;j<geometry_.numberOfSites();j++) {
				RealType value = jCoupling(i,j)*0.5;
				if (value==0) continue;
				value *= 0.5; // double counting i,j
				assert(i!=j);
				for (size_t orb2=0;orb2<mp_.orbitals;orb2++) {
					//if (orb2!=orb) continue; // testing only!!
					int sign = jTermSign(ket1,ket2,i,orb,j,orb2,basis);
					setSplusSminus(sparseRow,ket1,ket2,
							i,orb,j,orb2,value*sign,basis);
					setSplusSminus(sparseRow,ket1,ket2,
							j,orb2,i,orb,value*sign,basis);
				}
			}
		}

		int jTermSign(
				const WordType& ket1,
				const WordType& ket2,
				size_t i,
				size_t orb1,
				size_t j,
				size_t orb2,
				const BasisType &basis) const
		{
			if (i>j) return jTermSign(ket1,ket2,j,orb2,i,orb1,basis);
			int x = basis.doSign(ket1,ket2,i,orb1,j,orb2,SPIN_UP);
			x *= basis.doSign(ket1,ket2,i,orb1,j,orb2,SPIN_DOWN);
			return x;
		}

		void calcDiagonalElements(std::vector<RealType>& diag,
					  const BasisType &basis) const
		{
			size_t hilbert=basis.size();
			size_t nsite = geometry_.numberOfSites();

			// Calculate diagonal elements
			for (size_t ispace=0;ispace<hilbert;ispace++) {
				WordType ket1 = basis(ispace,SPIN_UP);
				WordType ket2 = basis(ispace,SPIN_DOWN);
				diag[ispace]=findS(nsite,ket1,ket2,ispace,basis);
			}
		}

		void calcDiagonalElements(VectorType &x,
					  const VectorType& y,
					  const BasisType* basis) const
		{
			size_t hilbert=basis->size();
			size_t nsite = geometry_.numberOfSites();

			for (size_t ispace=0;ispace<hilbert;ispace++) {
				WordType ket1 = basis->operator()(ispace,SPIN_UP);
				WordType ket2 = basis->operator()(ispace,SPIN_DOWN);
				x[ispace] += findS(nsite,ket1,ket2,ispace,*basis)*y[ispace];
			}
		}

		RealType findS(size_t nsite,WordType ket1,WordType ket2,size_t ispace,const BasisType& basis) const
		{
			RealType s = 0;
			for (size_t i=0;i<nsite;i++) {
				for (size_t orb=0;orb<mp_.orbitals;orb++) {

					// Hubbard term U0
					s += mp_.hubbardU[0] * basis.isThereAnElectronAt(ket1,ket2,
											 i,SPIN_UP,orb) * basis.isThereAnElectronAt(ket1,ket2,
																    i,SPIN_DOWN,orb);


					for (size_t orb2=orb+1;orb2<mp_.orbitals;orb2++) {
						// Hubbard term U1
						s += mp_.hubbardU[1] * nix(ket1,ket2,i,orb,basis) *
								nix(ket1,ket2,i,orb2,basis);

						// Diagonal U2 term
						s+= mp_.hubbardU[2]*
								szTerm(ket1,ket2,i,orb,basis)*
								szTerm(ket1,ket2,i,orb2,basis);
					}

					// JNN and JNNN diagonal part
					for (size_t j=0;j<nsite;j++) {
						for (size_t orb2=0;orb2<mp_.orbitals;orb2++) {
							RealType value = jCoupling(i,j);
							if (value==0) continue;
							s += value*0.5* // double counting i,j
									szTerm(ket1,ket2,i,orb,basis)*
									szTerm(ket1,ket2,j,orb2,basis);
						}
					}

					// Potential term
					s += mp_.potentialV[i+orb*nsite]*
							(basis.getN(ket1,i,SPIN_UP,orb) +
							 basis.getN(ket2,i,SPIN_DOWN,orb));

				}
			}
			return s;
		}

		size_t splusSminusNonZero(
						const WordType& ket1,
						const WordType& ket2,
						size_t i,
						size_t orb1,
						size_t j,
						size_t orb2,
						const BasisType &basis) const
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

		size_t u3TermNonZero(
				const WordType& ket1,
				const WordType& ket2,
				size_t i,
				size_t orb1,
				size_t orb2,
				const BasisType &basis) const
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

		size_t nix(
				const WordType& ket1,
				const WordType& ket2,
				size_t i,
				size_t orb,
				const BasisType &basis) const
		{
			size_t sum = 0;
			for (size_t spin=0;spin<2;spin++)
				sum += basis.isThereAnElectronAt(ket1,ket2,i,spin,orb);
			return sum;
		}

		RealType szTerm(
				const WordType& ket1,
				const WordType& ket2,
				size_t i,
				size_t orb,
				const BasisType &basis) const
		{
			RealType sz = basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb);
			sz -= basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb);
			return 0.5*sz;
		}
		
		RealType jCoupling(size_t i,size_t j) const
		{
			if (geometry_.terms()==1) return 0;
			return geometry_(i,0,j,0,TERM_J);
		}

		const ParametersModelType& mp_;
		const GeometryType& geometry_;
		BasisType basis_;
		
	}; // class FeBasedSc
	
} // namespace LanczosPlusPlus
#endif

