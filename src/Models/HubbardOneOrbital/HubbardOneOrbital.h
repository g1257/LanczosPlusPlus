
/*
*/

#ifndef HUBBARDLANCZOS_H
#define HUBBARDLANCZOS_H

#include "CrsMatrix.h"
#include "BasisHubbardLanczos.h"
#include "BitManip.h"
#include "TypeToString.h"
#include "SparseRow.h"
#include "ParametersModelHubbard.h"
#include "ProgramGlobals.h"

namespace LanczosPlusPlus {

	template<typename RealType_,typename GeometryType_>
	class HubbardOneOrbital {

		typedef PsimagLite::Matrix<RealType_> MatrixType;

	public:

		typedef ParametersModelHubbard<RealType_> ParametersModelType;
		typedef GeometryType_ GeometryType;
		typedef PsimagLite::CrsMatrix<RealType_> SparseMatrixType;
		typedef PsimagLite::SparseRow<SparseMatrixType> SparseRowType;
		typedef BasisHubbardLanczos<GeometryType> BasisType;
		typedef typename BasisType::WordType WordType;
		typedef RealType_ RealType;
		typedef std::vector<RealType> VectorType;

		static int const FERMION_SIGN = BasisType::FERMION_SIGN;

		HubbardOneOrbital(size_t nup,
		                  size_t ndown,
						  const ParametersModelType& mp,
						  const GeometryType& geometry)
		: mp_(mp),
		  geometry_(geometry),
		  basis_(geometry,nup,ndown),
		  hoppings_(geometry_.numberOfSites(),geometry_.numberOfSites())
		{
			size_t n = geometry_.numberOfSites();
			for (size_t i=0;i<n;i++)
				for (size_t j=0;j<n;j++)
					hoppings_(i,j) = geometry_(i,0,j,0,0);
		}

		size_t size() const { return basis_.size(); }

		size_t orbitals(size_t site) const
		{
			return 1;
		}

		void setupHamiltonian(SparseMatrixType &matrix) const
		{
			setupHamiltonian(matrix,basis_);
		}

		//! Gf. related functions below:
		void setupHamiltonian(SparseMatrixType &matrix,
		                      const BasisType &basis) const
		{
			size_t hilbert=basis.size();
			std::vector<RealType> diag(hilbert);
			calcDiagonalElements(diag,basis);
			
			size_t nsite = geometry_.numberOfSites();

			matrix.resize(hilbert,hilbert);
			// Calculate off-diagonal elements AND store matrix
			size_t nCounter=0;
			for (size_t ispace=0;ispace<hilbert;ispace++) {
				SparseRowType sparseRow;
				matrix.setRow(ispace,nCounter);
				WordType ket1 = basis(ispace,ProgramGlobals::SPIN_UP);
				WordType ket2 = basis(ispace,ProgramGlobals::SPIN_DOWN);
				// Save diagonal
				sparseRow.add(ispace,diag[ispace]);
				for (size_t i=0;i<nsite;i++) {
					setHoppingTerm(sparseRow,ket1,ket2,i,basis);
				}
				nCounter += sparseRow.finalize(matrix);
			}
			matrix.setRow(hilbert,nCounter);
		}

		bool hasNewParts(std::pair<size_t,size_t>& newParts,
		                 size_t what,
		                 size_t spin,
		                 const std::pair<size_t,size_t>& orbs) const
		{
			if (what==ProgramGlobals::OPERATOR_C || what==ProgramGlobals::OPERATOR_CDAGGER)
				return hasNewPartsCorCdagger(newParts,what,spin,orbs);
			std::string str(__FILE__);
			str += " " + ttos(__LINE__) +  "\n";
			str += std::string("hasNewParts: unsupported operator ");
			str += ProgramGlobals::id2Operator(what) + "\n";
			throw std::runtime_error(str.c_str());
		}

		const GeometryType& geometry() const { return geometry_; }

		const BasisType& basis() const { return basis_; }

		std::string name() const { return __FILE__; }

	private:

		bool hasNewPartsCorCdagger(std::pair<size_t,size_t>& newParts,
		                           size_t what,
		                           size_t spin,
		                           const std::pair<size_t,size_t>& orbs) const
		{
			int newPart1=basis_.electrons(ProgramGlobals::SPIN_UP);
			int newPart2=basis_.electrons(ProgramGlobals::SPIN_DOWN);
			int c = (what==ProgramGlobals::OPERATOR_C) ? -1 : 1;
			if (spin==ProgramGlobals::SPIN_UP) newPart1 += c;
			else newPart2 += c;

			if (newPart1<0 || newPart2<0) return false;
			size_t nsite = geometry_.numberOfSites();
			if (size_t(newPart1)>nsite || size_t(newPart2)>nsite) return false;
			if (newPart1==0 && newPart2==0) return false;
			newParts.first = size_t(newPart1);
			newParts.second = size_t(newPart2);
			return true;
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

					// Hubbard term U0
					s += mp_.hubbardU[i] *
							basis.isThereAnElectronAt(ket1,ket2,i,ProgramGlobals::SPIN_UP) *
							basis.isThereAnElectronAt(ket1,ket2,i,ProgramGlobals::SPIN_DOWN);

					// Potential term
					RealType tmp = mp_.potentialV[i];
					if (mp_.potentialT.size()>0) tmp += mp_.potentialT[i]*mp_.timeFactor;
					if (tmp!=0)
						s += tmp*
								(basis.getN(ket1,ket2,i,ProgramGlobals::SPIN_UP) +
								 basis.getN(ket1,ket2,i,ProgramGlobals::SPIN_DOWN));
				}
				diag[ispace]=s;
			}
		}

		void setHoppingTerm(SparseRowType &sparseRow,
		                    const WordType& ket1,
		                    const WordType& ket2,
		                    size_t i,
		                    const BasisType &basis) const
		{
			WordType s1i=(ket1 & BasisType::bitmask(i));
			if (s1i>0) s1i=1;
			WordType s2i=(ket2 & BasisType::bitmask(i));
			if (s2i>0) s2i=1;

			size_t nsite = geometry_.numberOfSites();

			// Hopping term
			for (size_t j=0;j<nsite;j++) {
				if (j<i) continue;
				RealType h = hoppings_(i,j);
				if (h==0) continue;
				WordType s1j= (ket1 & BasisType::bitmask(j));
				if (s1j>0) s1j=1;
				WordType s2j= (ket2 & BasisType::bitmask(j));
				if (s2j>0) s2j=1;

				if (s1i+s1j==1) {
					WordType bra1= ket1 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
					size_t temp = basis.perfectIndex(bra1,ket2);
					int extraSign = (s1i==1) ? FERMION_SIGN : 1;
					RealType cTemp = h*extraSign*basis_.doSign(ket1,ket2,i,j,ProgramGlobals::SPIN_UP);
					assert(temp<basis_.size());
					sparseRow.add(temp,cTemp);
				}

				if (s2i+s2j==1) {
					WordType bra2= ket2 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
					size_t temp = basis.perfectIndex(ket1,bra2);
					int extraSign = (s2i==1) ? FERMION_SIGN : 1;
					RealType cTemp = h*extraSign*basis_.doSign(ket1,ket2,i,j,ProgramGlobals::SPIN_DOWN);
					assert(temp<basis_.size());
					sparseRow.add(temp,cTemp);
				}
			}
		}

		const ParametersModelType& mp_;
		const GeometryType& geometry_;
		BasisType basis_;
		PsimagLite::Matrix<RealType> hoppings_;

	}; // class HubbardOneOrbital 
} // namespace LanczosPlusPlus
#endif

