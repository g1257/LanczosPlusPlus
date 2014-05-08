
/*
*/

#ifndef HUBBARDLANCZOS_H
#define HUBBARDLANCZOS_H

#include "BasisHubbardLanczos.h"
#include "BitManip.h"
#include "TypeToString.h"
#include "SparseRow.h"
#include "ParametersModelHubbard.h"
#include "ProgramGlobals.h"
#include "ModelBase.h"

namespace LanczosPlusPlus {

	template<typename RealType,typename GeometryType,typename InputType>
	class HubbardOneOrbital : public ModelBase<RealType,GeometryType,InputType> {

		typedef PsimagLite::Matrix<RealType> MatrixType;
		typedef ModelBase<RealType,GeometryType,InputType> BaseType;

		enum {TERM_HOPPING=0,TERM_NINJ=1,TERM_SUPER=2};

		enum {SPIN_UP = ProgramGlobals::SPIN_UP, SPIN_DOWN = ProgramGlobals::SPIN_DOWN};

	public:

		typedef ParametersModelHubbard<RealType,InputType> ParametersModelType;
		typedef BasisHubbardLanczos<GeometryType> BasisType;
		typedef typename BasisType::BaseType BasisBaseType;
                typedef typename BasisType::WordType WordType;
                typedef typename BaseType::SparseMatrixType SparseMatrixType;
                typedef typename BaseType::VectorType VectorType;
		typedef PsimagLite::SparseRow<SparseMatrixType> SparseRowType;

		static int const FERMION_SIGN = BasisType::FERMION_SIGN;

		HubbardOneOrbital(SizeType nup,
		                  SizeType ndown,
						  const ParametersModelType& mp,
						  const GeometryType& geometry)
		: mp_(mp),
		  geometry_(geometry),
		  basis_(geometry,nup,ndown),
		  hoppings_(geometry_.numberOfSites(),geometry_.numberOfSites())
		{
			SizeType n = geometry_.numberOfSites();
			for (SizeType i=0;i<n;i++)
				for (SizeType j=0;j<n;j++)
					hoppings_(i,j) = geometry_(i,0,j,0,TERM_HOPPING);
		}

		~HubbardOneOrbital()
		{
			BaseType::deleteGarbage(garbage_);
		}

		SizeType size() const { return basis_.size(); }

		SizeType orbitals(SizeType site) const
		{
			return 1;
		}

		void setupHamiltonian(SparseMatrixType& matrix) const
		{
			setupHamiltonian(matrix,basis_);
		}

		//! Gf. related functions below:
		void setupHamiltonian(SparseMatrixType& matrix,
		                      const BasisBaseType& basis) const
		{
			SizeType hilbert=basis.size();
			typename PsimagLite::Vector<RealType>::Type diag(hilbert);
			calcDiagonalElements(diag,basis);
			
			SizeType nsite = geometry_.numberOfSites();

			matrix.resize(hilbert,hilbert);
			// Calculate off-diagonal elements AND store matrix
			SizeType nCounter=0;
			for (SizeType ispace=0;ispace<hilbert;ispace++) {
				SparseRowType sparseRow;
				matrix.setRow(ispace,nCounter);
				WordType ket1 = basis(ispace,ProgramGlobals::SPIN_UP);
				WordType ket2 = basis(ispace,SPIN_DOWN);
				// Save diagonal
				sparseRow.add(ispace,diag[ispace]);
				for (SizeType i=0;i<nsite;i++) {
					setHoppingTerm(sparseRow,ket1,ket2,i,basis);
					setJTermOffDiagonal(sparseRow,ket1,ket2,i,basis);
				}

				nCounter += sparseRow.finalize(matrix);
			}
			matrix.setRow(hilbert,nCounter);
		}

		bool hasNewParts(std::pair<SizeType,SizeType>& newParts,
		                 SizeType what,
		                 SizeType spin,
		                 const std::pair<SizeType,SizeType>& orbs) const
		{
			if (what==ProgramGlobals::OPERATOR_C || what==ProgramGlobals::OPERATOR_CDAGGER)
				return hasNewPartsCorCdagger(newParts,what,spin,orbs);
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) +  "\n";
			str += PsimagLite::String("hasNewParts: unsupported operator ");
			str += ProgramGlobals::id2Operator(what) + "\n";
			throw std::runtime_error(str.c_str());
		}

		const GeometryType& geometry() const { return geometry_; }

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

		bool hasNewPartsCorCdagger(std::pair<SizeType,SizeType>& newParts,
		                           SizeType what,
		                           SizeType spin,
		                           const std::pair<SizeType,SizeType>& orbs) const
		{
			int newPart1=basis_.electrons(ProgramGlobals::SPIN_UP);
			int newPart2=basis_.electrons(SPIN_DOWN);
			int c = (what==ProgramGlobals::OPERATOR_C) ? -1 : 1;
			if (spin==ProgramGlobals::SPIN_UP) newPart1 += c;
			else newPart2 += c;

			if (newPart1<0 || newPart2<0) return false;
			SizeType nsite = geometry_.numberOfSites();
			if (SizeType(newPart1)>nsite || SizeType(newPart2)>nsite) return false;
			if (newPart1==0 && newPart2==0) return false;
			newParts.first = SizeType(newPart1);
			newParts.second = SizeType(newPart2);
			return true;
		}

		void calcDiagonalElements(typename PsimagLite::Vector<RealType>::Type& diag,
		                          const BasisBaseType& basis) const
		{
			SizeType hilbert=basis.size();
			SizeType nsite = geometry_.numberOfSites();
			SizeType orb = 0;

			// Calculate diagonal elements
			for (SizeType ispace=0;ispace<hilbert;ispace++) {
				WordType ket1 = basis(ispace,ProgramGlobals::SPIN_UP);
				WordType ket2 = basis(ispace,SPIN_DOWN);
				RealType s=0;
				for (SizeType i=0;i<nsite;i++) {

					// Hubbard term U0
					s += mp_.hubbardU[i] *
							basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb) *
							basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb);

					// SzSz
					for (SizeType j=0;j<nsite;j++) {
						RealType value = jCoupling(i,j);
						if (value==0) continue;
						s += value*0.5* // double counting i,j
						        szTerm(ket1,ket2,i,basis)*
						        szTerm(ket1,ket2,j,basis);
					}

					// Coulomb
					SizeType ne = (basis.getN(ket1,ket2,i,SPIN_UP,orb) +
					               basis.getN(ket1,ket2,i,SPIN_DOWN,orb));

					for (SizeType j=0;j<nsite;j++) {
						RealType value = coulombCoupling(i,j);
						if (value==0) continue;
						s += value * ne *
						        (basis.getN(ket1,ket2,j,SPIN_UP,orb) +
						         basis.getN(ket1,ket2,j,SPIN_DOWN,orb));
					}

					// Potential term
					RealType tmp = mp_.potentialV[i];
					if (mp_.potentialT.size()>0)
						tmp += mp_.potentialT[i]*mp_.timeFactor;
					if (tmp!=0) s += tmp * ne;
				}

				diag[ispace]=s;
			}
		}

		void setHoppingTerm(SparseRowType& sparseRow,
		                    const WordType& ket1,
		                    const WordType& ket2,
		                    SizeType i,
		                    const BasisBaseType& basis) const
		{
			WordType s1i=(ket1 & BasisType::bitmask(i));
			if (s1i>0) s1i=1;
			WordType s2i=(ket2 & BasisType::bitmask(i));
			if (s2i>0) s2i=1;

			SizeType nsite = geometry_.numberOfSites();
			SizeType orb = 0;

			// Hopping term
			for (SizeType j=0;j<nsite;j++) {
				if (j<i) continue;
				RealType h = hoppings_(i,j);
				if (h==0) continue;
				WordType s1j= (ket1 & BasisType::bitmask(j));
				if (s1j>0) s1j=1;
				WordType s2j= (ket2 & BasisType::bitmask(j));
				if (s2j>0) s2j=1;

				if (s1i+s1j==1) {
					WordType bra1= ket1 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
					SizeType temp = basis.perfectIndex(bra1,ket2);
					int extraSign = (s1i==1) ? FERMION_SIGN : 1;
					RealType cTemp = h*extraSign*
					                 basis.doSign(ket1,ket2,i,orb,j,orb,SPIN_UP);
					assert(temp<basis.size());
					sparseRow.add(temp,cTemp);
				}

				if (s2i+s2j==1) {
					WordType bra2= ket2 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
					SizeType temp = basis.perfectIndex(ket1,bra2);
					int extraSign = (s2i==1) ? FERMION_SIGN : 1;
					RealType cTemp = h*extraSign*basis.doSign(ket1,ket2,i,orb,j,orb,SPIN_DOWN);
					assert(temp<basis.size());
					sparseRow.add(temp,cTemp);
				}
			}
		}

		void setJTermOffDiagonal(SparseRowType& sparseRow,
		                         const WordType& ket1,
		                         const WordType& ket2,
		                         SizeType i,
		                         const BasisBaseType& basis) const
		{
			for (SizeType j=0;j<geometry_.numberOfSites();j++) {
				RealType value = jCoupling(i,j)*0.5;
				if (value==0) continue;
				value *= 0.5; // double counting i,j
				assert(i!=j);

				int sign = jTermSign(ket1,ket2,i,j,basis);
				setSplusSminus(sparseRow,ket1,ket2,i,j,value*sign,basis);
				setSplusSminus(sparseRow,ket1,ket2,j,i,value*sign,basis);
			}
		}

		void setSplusSminus(SparseRowType &sparseRow,
		                    const WordType& ket1,
		                    const WordType& ket2,
		                    SizeType i,
		                    SizeType j,
		                    RealType value,
		                    const BasisBaseType &basis) const
		{
			if (splusSminusNonZero(ket1,ket2,i,j,basis)==0) return;

			assert(i!=j);
			WordType bra1 = ket1 ^ (BasisType::bitmask(i)|BasisType::bitmask(j));
			WordType bra2 = ket2 ^ (BasisType::bitmask(i)|BasisType::bitmask(j));
			SizeType temp = basis.perfectIndex(bra1,bra2);
			sparseRow.add(temp,value);
		}

		SizeType splusSminusNonZero(const WordType& ket1,
		                            const WordType& ket2,
		                            SizeType i,
		                            SizeType j,
		                            const BasisBaseType &basis) const
		{
			SizeType orb = 0;
			if (basis.isThereAnElectronAt(ket1,ket2,j,SPIN_UP,orb)==0) return 0;
			if (basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb)==1) return 0;
			if (basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb)==0) return 0;
			if (basis.isThereAnElectronAt(ket1,ket2,j,SPIN_DOWN,orb)==1) return 0;
			return 1;
		}

		int jTermSign(const WordType& ket1,
		              const WordType& ket2,
		              SizeType i,
		              SizeType j,
		              const BasisBaseType& basis) const
		{
			if (i>j) return jTermSign(ket1,ket2,j,i,basis);
			SizeType orb = 0;
			int x = basis.doSign(ket1,ket2,i,orb,j,orb,SPIN_UP);
			x *= basis.doSign(ket1,ket2,i,orb,j,orb,SPIN_DOWN);
			return x;
		}


		RealType szTerm(const WordType& ket1,
		                const WordType& ket2,
		                SizeType i,
		                const BasisBaseType& basis) const
		{
			SizeType orb = 0;
			RealType sz = basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb);
			sz -= basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb);
			return 0.5*sz;
		}

		RealType jCoupling(SizeType i,SizeType j) const
		{
			if (geometry_.terms()<3) return 0;
			return geometry_(i,0,j,0,TERM_SUPER);
		}

		RealType coulombCoupling(SizeType i,SizeType j) const
		{
			if (geometry_.terms()<2) return 0;
			return geometry_(i,0,j,0,TERM_NINJ);
		}

		const ParametersModelType mp_;
		const GeometryType& geometry_;
		BasisType basis_;
		PsimagLite::Matrix<RealType> hoppings_;
		mutable typename PsimagLite::Vector<BasisType*>::Type garbage_;
	}; // class HubbardOneOrbital 
} // namespace LanczosPlusPlus
#endif

