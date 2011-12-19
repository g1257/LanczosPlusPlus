
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

#ifndef IMMM_HEADER_H
#define IMMM_HEADER_H

#include "CrsMatrix.h"
#include "BasisImmm.h"
#include "SparseRow.h"

namespace LanczosPlusPlus {

	template<typename RealType_,typename ParametersType,typename GeometryType>
	class Immm {

		typedef PsimagLite::Matrix<RealType_> MatrixType;

	public:

		typedef BasisImmm<GeometryType> BasisType;
		typedef typename BasisType::WordType WordType;
		typedef RealType_ RealType;
		typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;
		typedef PsimagLite::SparseRow<SparseMatrixType> SparseRowType;
		typedef std::vector<RealType> VectorType;

		enum {SPIN_UP=BasisType::SPIN_UP,SPIN_DOWN=BasisType::SPIN_DOWN};

		enum {DESTRUCTOR=BasisType::DESTRUCTOR,CONSTRUCTOR=BasisType::CONSTRUCTOR};

		static int const FERMION_SIGN = BasisType::FERMION_SIGN;

		Immm(size_t nup,size_t ndown,const ParametersType& mp,GeometryType& geometry)
		: mp_(mp),geometry_(geometry),basis_(geometry,nup,ndown)
		{}

		size_t size() const { return basis_.size(); }

		void setupHamiltonian(SparseMatrixType &matrix) const
		{
			setupHamiltonian(matrix,basis_);
		}

		bool hasNewParts(std::pair<size_t,size_t>& newParts,
		                 size_t type,
		                 size_t spin) const
		{
			int newPart1=basis_.electrons(SPIN_UP);
			int newPart2=basis_.electrons(SPIN_DOWN);
			int c = (type&1) ? -1 : 1;
			if (spin==SPIN_UP) newPart1 += c;
			else newPart2 += c;

			if (newPart1<0 || newPart2<0) return false;
			size_t nsite = geometry_.numberOfSites();
			if (size_t(newPart1)>nsite || size_t(newPart2)>nsite) return false;
			if (newPart1==0 && newPart2==0) return false;
			newParts.first = size_t(newPart1);
			newParts.second = size_t(newPart2);
			return true;
		}

		void getModifiedState(std::vector<RealType>& modifVector,
		                      const std::vector<RealType>& gsVector,
		                      const BasisType& basisNew,
		                      size_t type,
		                      size_t isite,
		                      size_t jsite,
		                      size_t spin) const
		{
			size_t what= (type&1) ?  DESTRUCTOR : CONSTRUCTOR;

			modifVector.resize(basisNew.size());
			for (size_t temp=0;temp<modifVector.size();temp++)
				modifVector[temp]=0.0;

			accModifiedState(modifVector,basisNew,gsVector,what,isite,spin,1);
			std::cerr<<"isite="<<isite<<" type="<<type;
			std::cerr<<" modif="<<(modifVector*modifVector)<<"\n";

			int isign= (type>1) ? -1 : 1;
			accModifiedState(modifVector,basisNew,gsVector,what,jsite,spin,isign);
			std::cerr<<"jsite="<<jsite<<" type="<<type;
			std::cerr<<" modif="<<(modifVector*modifVector)<<"\n";
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
			std::vector<RealType> diag(hilbert);
			calcDiagonalElements(diag,*basis);

			size_t nsite = geometry_.numberOfSites();

			// Calculate off-diagonal elements AND store matrix
			for (size_t ispace=0;ispace<hilbert;ispace++) {
				SparseRowType sparseRow;
				WordType ket1 = basis->operator()(ispace,SPIN_UP);
				WordType ket2 = basis->operator()(ispace,SPIN_DOWN);
				// Save diagonal
				sparseRow.add(ispace,diag[ispace]);
				for (size_t i=0;i<nsite;i++) {
					for (size_t orb=0;orb<basis->orbsPerSite(i);orb++) {
						setHoppingTerm(sparseRow,ket1,ket2,i,orb,*basis);
// 						if (orb==0) {
// 							setU2OffDiagonalTerm(sparseRow,ket1,ket2,
// 								i,orb,basis);
// 						}
// 						setU3Term(sparseRow,ket1,ket2,
// 								i,orb,1-orb,basis);
// 						setJTermOffDiagonal(sparseRow,ket1,ket2,
// 								i,orb,basis);
					}
				}
				//nCounter += sparseRow.finalize(matrix);
				x[ispace] += sparseRow.matrixVectorProduct(y);
			}
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
			matrix.resize(hilbert);

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
					for (size_t orb=0;orb<basis.orbsPerSite(i);orb++) {
						setHoppingTerm(sparseRow,ket1,ket2,i,orb,basis);
// 						if (orb==0) {
// 							setU2OffDiagonalTerm(sparseRow,ket1,ket2,
// 								i,orb,basis);
// 						}
// 						setU3Term(sparseRow,ket1,ket2,
// 								i,orb,1-orb,basis);
// 						setJTermOffDiagonal(sparseRow,ket1,ket2,
// 								i,orb,basis);
					}
				}
				nCounter += sparseRow.finalize(matrix);
			}
			matrix.setRow(hilbert,nCounter);
		}

	private:

		RealType hoppings(size_t i,size_t orb1,size_t j,size_t orb2) const
		{
			return geometry_(i,orb1,j,orb2,0);
		}

		void setHoppingTerm(
				SparseRowType &sparseRow,
				const WordType& ket1,
				const WordType& ket2,
				size_t i,
				size_t orb,
				const BasisType &basis) const
		{
			size_t ii = i*basis.orbs()+orb;
			WordType s1i=(ket1 & BasisType::bitmask(ii));
			if (s1i>0) s1i=1;
			WordType s2i=(ket2 & BasisType::bitmask(ii));
			if (s2i>0) s2i=1;

			size_t nsite = geometry_.numberOfSites();

			// Hopping term
			for (size_t j=0;j<nsite;j++) {
				if (j<i) continue;
				for (size_t orb2=0;orb2<basis.orbsPerSite(j);orb2++) {
					size_t jj = j*basis.orbs()+orb2;
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
						RealType cTemp = h*extraSign*basis_.doSign(ket1,ket2,i,orb,j,orb2,SPIN_UP);
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

// 		void setU2OffDiagonalTerm(
// 				SparseRowType &sparseRow,
// 				const WordType& ket1,
// 				const WordType& ket2,
// 				size_t i,
// 				size_t orb,
// 				const BasisType &basis) const
// 		{
// 			RealType val = FERMION_SIGN * mp_.hubbardU[2]*0.5;
// 			setSplusSminus(sparseRow,ket1,ket2,i,orb,i,1-orb,val,basis);
// 			setSplusSminus(sparseRow,ket1,ket2,i,1-orb,i,orb,val,basis);
// 		}
// 
// 		// N.B.: orb1!=orb2 here
// 		void setSplusSminus(
// 				SparseRowType &sparseRow,
// 				const WordType& ket1,
// 				const WordType& ket2,
// 				size_t i,
// 				size_t orb1,
// 				size_t j,
// 				size_t orb2,
// 				RealType value,
// 				const BasisType &basis) const
// 		{
// 			if (splusSminusNonZero(ket1,ket2,i,orb1,j,orb2,basis)==0) return;
// 
// 			size_t ii = i*ORBITALS + orb1;
// 			size_t jj = j*ORBITALS + orb2;
// 			assert(ii!=jj);
// 			WordType bra1 = ket1 ^ (BasisType::bitmask(ii)|BasisType::bitmask(jj));
// 			WordType bra2 = ket2 ^ (BasisType::bitmask(ii)|BasisType::bitmask(jj));
// 			size_t temp = basis.perfectIndex(bra1,bra2);
// 			sparseRow.add(temp,value);
// 		}
// 
// 		// N.B.: orb1!=orb2 here
// 		void setU3Term(
// 				SparseRowType &sparseRow,
// 				const WordType& ket1,
// 				const WordType& ket2,
// 				size_t i,
// 				size_t orb1,
// 				size_t orb2,
// 				const BasisType &basis) const
// 		{
// 			assert(orb1!=orb2);
// 			if (u3TermNonZero(ket1,ket2,i,orb1,orb2,basis)==0) return;
// 
// 			size_t ii = i*ORBITALS + orb1;
// 			size_t jj = i*ORBITALS + orb2;
// 			WordType bra1 = ket1 ^ (BasisType::bitmask(ii)|BasisType::bitmask(jj));
// 			WordType bra2 = ket2 ^ (BasisType::bitmask(ii)|BasisType::bitmask(jj));
// 			size_t temp = basis.perfectIndex(bra1,bra2);
// 			sparseRow.add(temp,FERMION_SIGN * mp_.hubbardU[3]);
// 		}

// 		void setJTermOffDiagonal(
// 				SparseRowType &sparseRow,
// 				const WordType& ket1,
// 				const WordType& ket2,
// 				size_t i,
// 				size_t orb,
// 				const BasisType &basis) const
// 		{
// 			for (size_t j=0;j<geometry_.numberOfSites();j++) {
// 				RealType value = jCoupling(i,j)*0.5;
// 				if (value==0) continue;
// 				value *= 0.5; // double counting i,j
// 				assert(i!=j);
// 				for (size_t orb2=0;orb2<ORBITALS;orb2++) {
// 					//if (orb2!=orb) continue; // testing only!!
// 					int sign = jTermSign(ket1,ket2,i,orb,j,orb2,basis);
// 					setSplusSminus(sparseRow,ket1,ket2,
// 							i,orb,j,orb2,value*sign,basis);
// 					setSplusSminus(sparseRow,ket1,ket2,
// 							j,orb2,i,orb,value*sign,basis);
// 				}
// 			}
// 
// 		}

// 		int jTermSign(
// 				const WordType& ket1,
// 				const WordType& ket2,
// 				size_t i,
// 				size_t orb1,
// 				size_t j,
// 				size_t orb2,
// 				const BasisType &basis) const
// 		{
// 			if (i>j) return jTermSign(ket1,ket2,j,orb2,i,orb1,basis);
// 			int x = basis.doSign(ket1,ket2,i,orb1,j,orb2,SPIN_UP);
// 			x *= basis.doSign(ket1,ket2,i,orb1,j,orb2,SPIN_DOWN);
// 			return x;
// 		}

		void calcDiagonalElements(std::vector<RealType>& diag,
		                          const BasisType &basis) const
		{
			size_t hilbert=basis.size();
			size_t nsite = geometry_.numberOfSites();

			// Calculate diagonal elements
			for (size_t ispace=0;ispace<hilbert;ispace++) {
				WordType ket1 = basis(ispace,SPIN_UP);
				WordType ket2 = basis(ispace,SPIN_DOWN);
				RealType s=0;
				for (size_t i=0;i<nsite;i++) {
					for (size_t orb=0;orb<basis.orbsPerSite(i);orb++) {

						// Hubbard term U0
						s += mp_.hubbardU[0] * basis.isThereAnElectronAt(ket1,ket2,
								i,SPIN_UP,orb) * basis.isThereAnElectronAt(ket1,ket2,
								i,SPIN_DOWN,orb);
						
						// Hubbard term U1
						if (orb==0) {
							s += mp_.hubbardU[1] * nix(ket1,ket2,i,orb,basis) *
									nix(ket1,ket2,i,1-orb,basis);
						}

						// Diagonal U2 term
// 						if (orb==0 && mp_.hubbardU[2]!=0) {
// 							s+= mp_.hubbardU[2]*
// 								szTerm(ket1,ket2,i,orb,basis)*
// 								szTerm(ket1,ket2,i,1-orb,basis);
// 						}

// 						// JNN and JNNN diagonal part
// 						for (size_t j=0;j<nsite;j++) {
// 							for (size_t orb2=0;orb2<ORBITALS;orb2++) {
// 								RealType value = jCoupling(i,j);
// 								if (value==0) continue;
// 								s += value*0.5* // double counting i,j
// 									szTerm(ket1,ket2,i,orb,basis)*
// 									szTerm(ket1,ket2,j,orb2,basis);
// 							}
// 						}

						// Potential term
						if (mp_.potentialV[i]!=0)
							s += mp_.potentialV[i]*
								(basis.getN(ispace,SPIN_UP,orb) +
								basis.getN(ispace,SPIN_DOWN,orb));
					}
				}
				diag[ispace]=s;
			}
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

// 		RealType szTerm(
// 				const WordType& ket1,
// 				const WordType& ket2,
// 				size_t i,
// 				size_t orb,
// 				const BasisType &basis) const
// 		{
// 			RealType sz = basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb);
// 			sz -= basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb);
// 			return 0.5*sz;
// 		}
		
// 		RealType jCoupling(size_t i,size_t j) const
// 		{
// 			if (geometry_.terms()==1) return 0;
// 			return geometry_(i,0,j,0,TERM_J);
// 		}

		//! Gf Related function:
		void accModifiedState(std::vector<RealType> &z,
		                      const BasisType& newBasis,
		                      const std::vector<RealType> &gsVector,
		                      size_t what,
		                      size_t site,
		                      size_t spin,
		                      int isign) const
		{
			for (size_t orb=0;orb<newBasis.orbsPerSite(site);orb++)
				accModifiedState(z,newBasis,gsVector,what,site,spin,orb,isign);
		}

		//! Gf Related function:
		void accModifiedState(std::vector<RealType> &z,
		                      const BasisType& newBasis,
		                      const std::vector<RealType> &gsVector,
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
				int mysign = basis_.doSignGf(ket1,ket2,site,spin,orb);
				z[temp] += isign*mysign*gsVector[ispace];
			}
		}

		const ParametersType& mp_;
		const GeometryType& geometry_;
		BasisType basis_;
	}; // class Immm
	
} // namespace LanczosPlusPlus
#endif

