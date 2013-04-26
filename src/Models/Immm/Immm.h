
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
#include "SparseRowCached.h"
#include "ParametersImmm.h"

namespace LanczosPlusPlus {

template<typename RealType_,typename GeometryType_>
	class Immm {

		typedef PsimagLite::Matrix<RealType_> MatrixType;

	public:

		typedef GeometryType_ GeometryType;
		typedef ParametersImmm<RealType_> ParametersModelType;
		typedef BasisImmm<GeometryType> BasisType;
		typedef typename BasisType::WordType WordType;
		typedef RealType_ RealType;
		typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;
		typedef PsimagLite::SparseRowCached<SparseMatrixType> SparseRowType;
		typedef std::vector<RealType> VectorType;

		static int const FERMION_SIGN = BasisType::FERMION_SIGN;

		Immm(size_t nup,size_t ndown,const ParametersModelType& mp,const GeometryType& geometry)
		: mp_(mp),geometry_(geometry),basis_(geometry,nup,ndown)
		{}

		size_t size() const { return basis_.size(); }

		size_t orbitals(size_t site) const
		{
			return basis_.orbsPerSite(site);
		}

		void setupHamiltonian(SparseMatrixType &matrix) const
		{
			if (!mp_.useReflectionSymmetry) {
				setupHamiltonian(matrix,basis_);
				return;
			}
//			SparseMatrixType matrix2;
//			setupHamiltonian(matrix2,basis_);
//			ReflectionSymmetryType rs(basis_,geometry_);
//			rs.transform(matrix,matrix2);
		}

		bool hasNewParts(std::pair<size_t,size_t>& newParts,
						 size_t what,
						 size_t spin,
						 const std::pair<size_t,size_t>& orbs) const
		{
			return basis_.hasNewParts(newParts,what,spin,orbs);
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
			size_t cacheSize = hilbert/10;
			if (cacheSize<100) cacheSize=100;
			SparseRowType sparseRow(cacheSize);
			for (size_t ispace=0;ispace<hilbert;ispace++) {
				WordType ket1 = basis->operator()(ispace,ProgramGlobals::SPIN_UP);
				WordType ket2 = basis->operator()(ispace,ProgramGlobals::SPIN_DOWN);
				// Save diagonal
				sparseRow.add(ispace,diag[ispace]);
				for (size_t i=0;i<nsite;i++) {
					for (size_t orb=0;orb<basis->orbsPerSite(i);orb++) {
						setHoppingTerm(sparseRow,ket1,ket2,ispace,i,orb,*basis);
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

		const BasisType& basis() const { return basis_; }

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
			SparseRowType sparseRow(100);
			for (size_t ispace=0;ispace<hilbert;ispace++) {
				matrix.setRow(ispace,nCounter);
				WordType ket1 = basis(ispace,ProgramGlobals::SPIN_UP);
				WordType ket2 = basis(ispace,ProgramGlobals::SPIN_DOWN);
				// Save diagonal
				sparseRow.add(ispace,diag[ispace]);
				for (size_t i=0;i<nsite;i++) {
					for (size_t orb=0;orb<basis.orbsPerSite(i);orb++) {
						setHoppingTerm(sparseRow,ket1,ket2,ispace,i,orb,basis);
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

		std::string name() const { return __FILE__; }

	private:

		RealType hoppings(size_t i,size_t orb1,size_t j,size_t orb2) const
		{
			return geometry_(i,orb1,j,orb2,0);
		}

		RealType Upd(size_t i,size_t j) const
		{
			return geometry_(i,0,j,0,1);
		}

		void setHoppingTerm(SparseRowType &sparseRow,
		                    const WordType& ket1,
		                    const WordType& ket2,
		                    size_t ispace,
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
			for (size_t j=i;j<nsite;j++) {
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
						size_t temp = basis.perfectIndex(bra1,ispace,ProgramGlobals::SPIN_UP);
						int extraSign = (s1i==1) ? FERMION_SIGN : 1;
						RealType cTemp = h*extraSign*basis.doSign(ket1,ket2,i,orb,j,orb2,ProgramGlobals::SPIN_UP);
						sparseRow.add(temp,cTemp);
					}

					if (s2i+s2j==1) {
						WordType bra2= ket2 ^(BasisType::bitmask(ii)|BasisType::bitmask(jj));
						size_t temp = basis.perfectIndex(bra2,ispace,ProgramGlobals::SPIN_DOWN);
						int extraSign = (s2i==1) ? FERMION_SIGN : 1;
						RealType cTemp = h*extraSign*basis.doSign(
							ket1,ket2,i,orb,j,orb2,ProgramGlobals::SPIN_DOWN);
						sparseRow.add(temp,cTemp);
					}
				}
			}
		}

		void calcDiagonalElements(std::vector<RealType>& diag,
		                          const BasisType &basis) const
		{
			size_t hilbert=basis.size();
			size_t nsite = geometry_.numberOfSites();

			// Calculate diagonal elements
			for (size_t ispace=0;ispace<hilbert;ispace++) {
				WordType ket1 = basis(ispace,ProgramGlobals::SPIN_UP);
				WordType ket2 = basis(ispace,ProgramGlobals::SPIN_DOWN);
				RealType s=0;
				for (size_t i=0;i<nsite;i++) {
					for (size_t orb=0;orb<basis.orbsPerSite(i);orb++) {

						size_t totalCharge = basis.getN(ket1,i,ProgramGlobals::SPIN_UP,orb) +
								basis.getN(ket2,i,ProgramGlobals::SPIN_DOWN,orb);

						// Hubbard term U0
						s += mp_.hubbardU[i] * (1.0-basis.isThereAnElectronAt(ket1,ket2,
								i,ProgramGlobals::SPIN_UP,orb)) * (1.0-basis.isThereAnElectronAt(ket1,ket2,
								i,ProgramGlobals::SPIN_DOWN,orb));

						// Potential term
						s += mp_.potentialV[i]*totalCharge;

						// Upd n_O n_Cu
						if (basis.orbsPerSite(i)==1) continue;
						// i is an Oxygen site now
						for (size_t j=0;j<nsite;j++) {
							if (basis.orbsPerSite(j)==2) continue;
							// j is a Copper site now
							size_t totalCharge2 = basis.getN(ket1,j,ProgramGlobals::SPIN_UP,0) +
									basis.getN(ket2,j,ProgramGlobals::SPIN_DOWN,0);
							s += (2.0-totalCharge) * (2.0-totalCharge2) * Upd(i,j);
						}
					}
				}
				diag[ispace]=s;
			}
		}

		//! Gf Related function:
		template<typename SomeVectorType>
		void accModifiedState(SomeVectorType& z,
							  size_t operatorLabel,
		                      const BasisType& newBasis,
							  const SomeVectorType& gsVector,
		                      size_t site,
		                      size_t spin,
		                      int isign) const
		{
			for (size_t orb=0;orb<newBasis.orbsPerSite(site);orb++)
				accModifiedState(z,operatorLabel,newBasis,gsVector,site,spin,orb,isign);
		}

		const ParametersModelType& mp_;
		const GeometryType& geometry_;
		BasisType basis_;
	}; // class Immm

} // namespace LanczosPlusPlus
#endif

