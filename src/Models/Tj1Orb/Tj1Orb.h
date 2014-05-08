
/*
*/

#ifndef TJ_1ORB_H
#define TJ_1ORB_H

#include "CrsMatrix.h"
#include "BasisTj1OrbLanczos.h"
#include "BitManip.h"
#include "TypeToString.h"
#include "SparseRow.h"
#include "ParametersTj1Orb.h"
#include "ModelBase.h"

namespace LanczosPlusPlus {

template<typename RealType,typename GeometryType,typename InputType>
class Tj1Orb  : public ModelBase<RealType,GeometryType,InputType> {

	typedef PsimagLite::Matrix<RealType> MatrixType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef ModelBase<RealType,GeometryType,InputType> BaseType;

public:

	typedef ParametersTj1Orb<RealType,InputType> ParametersModelType;
	typedef BasisTj1OrbLanczos<GeometryType> BasisType;
	typedef typename BasisType::WordType WordType;
	typedef typename BaseType::SparseMatrixType SparseMatrixType;
	typedef typename BaseType::VectorType VectorType;
	typedef PsimagLite::SparseRow<SparseMatrixType> SparseRowType;

	static int const FERMION_SIGN = BasisType::FERMION_SIGN;

	Tj1Orb(SizeType nup,
	       SizeType ndown,
	       InputType& io,
	       const GeometryType& geometry)
	    : mp_(io),
	      geometry_(geometry),
	      basis_(geometry,nup,ndown),
	      hoppings_(geometry_.numberOfSites(),geometry_.numberOfSites()),
	      j_(geometry_.numberOfSites(),geometry_.numberOfSites()),
	      w_(geometry_.numberOfSites(),geometry_.numberOfSites())
	{
		SizeType n = geometry_.numberOfSites();
		bool h = (geometry_.terms() == 1);

		for (SizeType i=0;i<n;i++) {
			for (SizeType j=0;j<n;j++) {
				hoppings_(i,j) = (h) ? 0 : geometry_(i,0,j,0,0);
				j_(i,j) = (h) ? geometry_(i,0,j,0,0) : geometry_(i,0,j,0,1);
				w_(i,j) = (h) ? 0 : geometry_(i,0,j,0,2);
			}
		}
	}

	SizeType size() const { return basis_.size(); }

	SizeType orbitals(SizeType site) const
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
		SizeType hilbert=basis.size();
		typename PsimagLite::Vector<RealType>::Type diag(hilbert,0.0);
		calcDiagonalElements(diag,basis);

		SizeType nsite = geometry_.numberOfSites();

		matrix.resize(hilbert,hilbert);
		// Calculate off-diagonal elements AND store matrix
		SizeType nCounter=0;
		for (SizeType ispace=0;ispace<hilbert;ispace++) {
			SparseRowType sparseRow;
			matrix.setRow(ispace,nCounter);
			WordType ket1 = basis(ispace,ProgramGlobals::SPIN_UP);
			WordType ket2 = basis(ispace,ProgramGlobals::SPIN_DOWN);
			//				std::cout<<"ket1="<<ket1<<" ket2="<<ket2<<"\n";
			// Save diagonal
			sparseRow.add(ispace,diag[ispace]);
			for (SizeType i=0;i<nsite;i++) {
				setHoppingTerm(sparseRow,ket1,ket2,i,basis);
				setSplusSminus(sparseRow,ket1,ket2,i,basis);
			}
			nCounter += sparseRow.finalize(matrix);
		}
		matrix.setRow(hilbert,nCounter);
		matrix.checkValidity();
		assert(isHermitian(matrix));
	}

	bool hasNewParts(std::pair<SizeType,SizeType>& newParts,
	                 SizeType what,
	                 SizeType spin,
	                 const PairType& orbs) const
	{
		if (what==ProgramGlobals::OPERATOR_C || what==ProgramGlobals::OPERATOR_CDAGGER)
			return hasNewPartsCorCdagger(newParts,what,spin,orbs);

		if (what == ProgramGlobals::OPERATOR_SPLUS ||
		    what == ProgramGlobals::OPERATOR_SMINUS)
			return hasNewPartsSplusOrMinus(newParts,what,spin,orbs);

		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) +  "\n";
		str += PsimagLite::String("hasNewParts: unsupported operator ");
		str += ProgramGlobals::id2Operator(what) + "\n";
		throw std::runtime_error(str.c_str());
	}

	const GeometryType& geometry() const { return geometry_; }

	const BasisType& basis() const { return basis_; }

	void printBasis(std::ostream& os) const
	{
		os<<basis_;
	}

	PsimagLite::String name() const { return __FILE__; }

private:

	bool hasNewPartsCorCdagger(std::pair<SizeType,SizeType>& newParts,
	                           SizeType what,
	                           SizeType spin,
	                           const PairType& orbs) const
	{
		int newPart1=basis_.electrons(ProgramGlobals::SPIN_UP);
		int newPart2=basis_.electrons(ProgramGlobals::SPIN_DOWN);
		int c = (what==ProgramGlobals::OPERATOR_CDAGGER) ? 1 : -1;
		if (spin==ProgramGlobals::SPIN_UP) newPart1 += c;
		else newPart2 += c;

		if (newPart1<0 || newPart2<0) return false;
		SizeType nsite = geometry_.numberOfSites();
		if (SizeType(newPart1)>nsite || SizeType(newPart2)>nsite) return false;
		if (newPart1==0 && newPart2==0) return false;
		if (SizeType(newPart1+newPart2)>nsite) return false; // no double occupancy
		newParts.first = SizeType(newPart1);
		newParts.second = SizeType(newPart2);
		return true;
	}

	bool hasNewPartsSplusOrMinus(std::pair<SizeType,SizeType>& newParts,
	                             SizeType what,
	                             SizeType spin,
	                             const PairType& orbs) const
	{
		int newPart1=basis_.electrons(ProgramGlobals::SPIN_UP);
		int newPart2=basis_.electrons(ProgramGlobals::SPIN_DOWN);
		int c = (what==ProgramGlobals::OPERATOR_SPLUS) ? 1 : -1;
		if (spin==ProgramGlobals::SPIN_UP) {
			newPart1 += c;
			newPart2 -= c;
		} else {
			newPart2 += c;
			newPart1 -= c;
		}

		if (newPart1<0 || newPart2<0) return false;

		SizeType nsite = geometry_.numberOfSites();
		if (SizeType(newPart1)>nsite || SizeType(newPart2)>nsite) return false;
		if (newPart1==0 && newPart2==0) return false;
		if (SizeType(newPart1+newPart2)>nsite) return false; // no double occupancy
		newParts.first = SizeType(newPart1);
		newParts.second = SizeType(newPart2);
		return true;
	}

	void calcDiagonalElements(typename PsimagLite::Vector<RealType>::Type& diag,
	                          const BasisType &basis) const
	{
		SizeType hilbert=basis.size();
		SizeType nsite = geometry_.numberOfSites();

		// Calculate diagonal elements
		for (SizeType ispace=0;ispace<hilbert;ispace++) {
			WordType ket1 = basis(ispace,ProgramGlobals::SPIN_UP);
			WordType ket2 = basis(ispace,ProgramGlobals::SPIN_DOWN);
			RealType s=0;
			for (SizeType i=0;i<nsite;i++) {

				int niup = basis.isThereAnElectronAt(ket1,ket2,i,ProgramGlobals::SPIN_UP);

				int nidown = basis.isThereAnElectronAt(ket1,ket2,i,ProgramGlobals::SPIN_DOWN);

				if (i < mp_.potentialV.size()) {
					s += mp_.potentialV[i]*niup;
					s += mp_.potentialV[i]*nidown;
				}

				for (SizeType j=i+1;j<nsite;j++) {

					int njup = basis.isThereAnElectronAt(ket1,ket2,j,ProgramGlobals::SPIN_UP);
					int njdown = basis.isThereAnElectronAt(ket1,ket2,j,ProgramGlobals::SPIN_DOWN);

					// Sz Sz term:
					s += (niup-nidown) * (njup - njdown)  * j_(i,j)*0.25;

					// ni nj term
					s+= (niup+nidown) * (njup + njdown) * w_(i,j);
				}
			}
			diag[ispace]=s;
		}
	}

	void setHoppingTerm(SparseRowType &sparseRow,
	                    const WordType& ket1,
	                    const WordType& ket2,
	                    SizeType i,
	                    const BasisType &basis) const
	{
		WordType s1i=(ket1 & BasisType::bitmask(i));
		if (s1i>0) s1i=1;
		WordType s2i=(ket2 & BasisType::bitmask(i));
		if (s2i>0) s2i=1;

		SizeType nsite = geometry_.numberOfSites();

		// Hopping term
		for (SizeType j=0;j<nsite;j++) {
			if (j<i) continue;
			RealType h = hoppings_(i,j);
			if (h==0) continue;
			WordType s1j= (ket1 & BasisType::bitmask(j));
			if (s1j>0) s1j=1;
			WordType s2j= (ket2 & BasisType::bitmask(j));
			if (s2j>0) s2j=1;

			if (s1i+s1j==1 && !(s1j==0 && s2j>0) && !(s1j>0 && s2i>0)) {
				WordType bra1= ket1 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
				SizeType temp = basis.perfectIndex(bra1,ket2);
				int extraSign = (s1i==1) ? FERMION_SIGN : 1;
				RealType cTemp = h*extraSign*basis_.doSign(ket1,ket2,i,j,ProgramGlobals::SPIN_UP);
				sparseRow.add(temp,cTemp);
			}

			if (s2i+s2j==1 && !(s2j==0 && s1j>0) && !(s2j>0 && s1i>0)) {
				WordType bra2= ket2 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
				SizeType temp = basis.perfectIndex(ket1,bra2);
				int extraSign = (s2i==1) ? FERMION_SIGN : 1;
				RealType cTemp = h*extraSign*
				                 basis_.doSign(ket1,ket2,i,j,ProgramGlobals::SPIN_DOWN);
				sparseRow.add(temp,cTemp);
			}
		}
	}

	void setSplusSminus(SparseRowType &sparseRow,
	                    const WordType& ket1,
	                    const WordType& ket2,
	                    SizeType i,
	                    const BasisType &basis) const
	{
		WordType s1i=(ket1 & BasisType::bitmask(i));
		if (s1i>0) s1i=1;
		WordType s2i=(ket2 & BasisType::bitmask(i));
		if (s2i>0) s2i=1;

		SizeType nsite = geometry_.numberOfSites();

		// Hopping term
		for (SizeType j=0;j<nsite;j++) {
			if (j<i) continue;
			RealType h = j_(i,j)*0.5;
			if (h==0) continue;
			WordType s1j= (ket1 & BasisType::bitmask(j));
			if (s1j>0) s1j=1;
			WordType s2j= (ket2 & BasisType::bitmask(j));
			if (s2j>0) s2j=1;

			if (s1i==1 && s1j==0 && s2i==0 && s2j==1) {
				WordType bra1= ket1 ^ BasisType::bitmask(i);
				bra1 |= BasisType::bitmask(j);
				WordType bra2= ket2 | BasisType::bitmask(i);
				bra2 ^= BasisType::bitmask(j);
				SizeType temp = basis.perfectIndex(bra1,bra2);
				sparseRow.add(temp,h*signSplusSminus(i,j,bra1,bra2));
			}

			if (s1i==0 && s1j==1 && s2i==1 && s2j==0) {
				WordType bra1= ket1 | BasisType::bitmask(i);
				bra1 ^= BasisType::bitmask(j);
				WordType bra2= ket2 ^ BasisType::bitmask(i);
				bra2 |= BasisType::bitmask(j);
				SizeType temp = basis.perfectIndex(bra1,bra2);
				sparseRow.add(temp,h*signSplusSminus(i,j,bra1,bra2));
			}
		}
	}

	int signSplusSminus(SizeType i,
	                    SizeType j,
	                    const WordType& bra1,
	                    const WordType& bra2) const
	{
		SizeType n = geometry_.numberOfSites();
		int s = 1;
		if (j>0) s *= parityFrom(0,j-1,bra2);
		if (i>0) s *= parityFrom(0,i-1,bra2);
		if (i<n-1) s *= parityFrom(i+1,n-1,bra1);
		if (j<n-1) s *= parityFrom(j+1,n-1,bra1);
		return s;
	}

	// from i to j including i and j
	// assumes i<=j
	int parityFrom(SizeType i,SizeType j,const WordType& ket) const
	{
		if (i==j) return (BasisType::bitmask(j) & ket) ? -1 : 1;
		assert(i<j);
		//j>i>=0 now
		WordType mask = ket;
		mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
		int s=(PsimagLite::BitManip::count(mask) & 1) ? -1 : 1;
		// Is there something of this species at i?
		if (BasisType::bitmask(i) & ket) s = -s;
		// Is there something of this species at j?
		if (BasisType::bitmask(j) & ket) s = -s;
		return s;
	}

	const ParametersModelType mp_;
	const GeometryType& geometry_;
	BasisType basis_;
	PsimagLite::Matrix<RealType> hoppings_;
	PsimagLite::Matrix<RealType> j_;
	PsimagLite::Matrix<RealType> w_;

}; // class Tj1Orb
} // namespace LanczosPlusPlus
#endif

