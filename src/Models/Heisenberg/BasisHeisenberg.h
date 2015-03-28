/*
*/

#ifndef LANCZOS_BASIS_HEISENBERG_H
#define LANCZOS_BASIS_HEISENBERG_H

#include "BitManip.h"
#include "ProgramGlobals.h"
#include "../../Engine/BasisBase.h"

namespace LanczosPlusPlus {

template<typename GeometryType>
class BasisHeisenberg : public BasisBase<GeometryType> {

	typedef ProgramGlobals::PairIntType PairIntType;

public:

	typedef BasisBase<GeometryType> BaseType;
	typedef typename BaseType::WordType WordType;
	typedef typename BaseType::VectorWordType VectorWordType;
	static VectorWordType bitmask_;

	enum {OPERATOR_NIL=ProgramGlobals::OPERATOR_NIL,
	      OPERATOR_SZ=ProgramGlobals::OPERATOR_SZ};

	BasisHeisenberg(const GeometryType& geometry,
	                SizeType twiceS,
	                SizeType szPlusConst)
	    : geometry_(geometry),
	      twiceS_(twiceS),
	      szPlusConst_(szPlusConst),
	      bits_(0)
	{
		SizeType sites = geometry_.numberOfSites();
		assert(bitmask_.size()==0 || bitmask_.size()== sites);
		if (bitmask_.size()==0) doBitmask();
		WordType searchTotal = 1;
		assert(twiceS > 0);
		bits_ = 1 + static_cast<SizeType>(logBase2(twiceS + 1));
		if (twiceS & 1) bits_--;
		searchTotal <<= (bits_ * sites);

		WordType mask = getMask();

		for (WordType lui = 0; lui < searchTotal; ++lui) {
			if (mOf(lui,mask) != szPlusConst) continue;
			data_.push_back(lui);
		}

		std::cerr<<(*this);
	}

	static const WordType& bitmask(SizeType i)
	{
		return bitmask_[i];
	}

	SizeType size() const { return data_.size(); }

	SizeType dofs() const { return twiceS_ + 1; }

	SizeType perfectIndex(const VectorWordType& kets) const
	{
		throw PsimagLite::RuntimeError("BasisHeisenberg::perfectIndex kets\n");
	}

	SizeType perfectIndex(WordType ket,WordType) const
	{
		for (SizeType i = 0; i < data_.size(); ++i) {
			if (ket == data_[i]) return i;
		}

		throw PsimagLite::RuntimeError("perfectIndex: no index found\n");
	}

	WordType operator()(SizeType i, SizeType) const
	{
		return data_[i];
	}

	SizeType isThereAnElectronAt(WordType ket1,
	                             WordType ket2,
	                             SizeType site,
	                             SizeType spin,
	                             SizeType) const
	{
		throw PsimagLite::RuntimeError("BasisHeisenberg::isThereAnElectronAt\n");
	}

	SizeType getN(WordType ket1,
	              WordType,
	              SizeType site,
	              SizeType,
	              SizeType) const
	{
		WordType mask = getMask();
		ket1 >>= (bits_ * site);
		return (ket1 & mask);
	}

	int doSignGf(WordType a, WordType b,SizeType ind,SizeType sector,SizeType) const
	{
		throw PsimagLite::RuntimeError("BasisHeisenberg::doSignGf\n");
	}

	int doSign(WordType ket1,
	           WordType ket2,
	           SizeType i,
	           SizeType,
	           SizeType j,
	           SizeType,
	           SizeType spin) const
	{
		throw PsimagLite::RuntimeError("BasisHeisenberg::doSign\n");
	}

	PairIntType getBraIndex(WordType ket1,
	                        WordType ket2,
	                        SizeType operatorLabel,
	                        SizeType site,
	                        SizeType spin,
	                        SizeType orb) const
	{
		if (operatorLabel == ProgramGlobals::OPERATOR_SPLUS) {
			throw PsimagLite::RuntimeError("BasisHeisenberg::getBraIndex SPLUS\n");

		} else if (operatorLabel == ProgramGlobals::OPERATOR_SMINUS) {
			throw PsimagLite::RuntimeError("BasisHeisenberg::getBraIndex SMINUS\n");
		}

		return getBraIndex_(ket1,ket2,operatorLabel,site,spin,orb);
	}

	SizeType orbsPerSite(SizeType) const { return 1; }

	SizeType orbs() const { return 1; }

	SizeType perfectIndex(WordType,SizeType,SizeType) const
	{
		throw PsimagLite::RuntimeError("BasisHeisenberg::perfectIndex\n");
	}

	void print(std::ostream& os) const
	{
		SizeType hilbert = 1;
		hilbert <<= geometry_.numberOfSites();
		ProgramGlobals::printBasisVector(os,hilbert,data_);
	}

	WordType getBra(WordType ket,
	                SizeType site1,
	                SizeType val1,
	                SizeType site2,
	                SizeType val2) const
	{
		WordType bra = ket;
		WordType mask1 = getMask();
		WordType mask2 = mask1;
		mask1 <<= (site1*bits_);
		bra &= (~mask1);

		mask2 <<= (site2*bits_);
		bra &= (~mask2);

		mask1 = val1;
		mask1 <<= (site1*bits_);
		bra |= mask1;

		mask2 = val2;
		mask2 <<= (site2*bits_);
		bra |= mask2;

		return bra;
	}

	template<typename GeometryType2>
	friend std::ostream& operator<<(std::ostream& os,
	                                const BasisHeisenberg<GeometryType2>& basis);

private:

	WordType getMask() const
	{
		SizeType mask = 1;
		for (SizeType i = 0; i < bits_; ++i)
			mask |= bitmask_[i];
		return mask;
	}

	int mOf(WordType lui, WordType mask) const
	{
		SizeType m = 0;
		bool allowed = true;
		while (lui != 0) {
			WordType tmp = (lui & mask);
			if (!mOfIsAllowed(tmp)) {
				allowed = false;
				break;
			}

			m += tmp;
			lui >>= bits_;
		}

		return (allowed) ? m : -1;
	}

	bool mOfIsAllowed(WordType val) const
	{
		if (twiceS_ & 1) return true;

		return (val <= twiceS_);
	}

	PairIntType getBraIndex_(const WordType& ket1,
	                         const WordType& ket2,
	                         SizeType operatorLabel,
	                         SizeType site,
	                         SizeType spin,
	                         SizeType) const
	{
		throw PsimagLite::RuntimeError("BasisHeisenberg::getBraIndex_ \n");
	}

	int getBra(WordType& bra,
	           SizeType operatorLabel,
	           const WordType& ket1,
	           const WordType& ket2,
	           //		            SizeType what,
	           SizeType site,
	           SizeType spin) const
	{
		switch(operatorLabel) {
		case ProgramGlobals::OPERATOR_SZ:
			return getBraSzOrN(bra,ket1,ket2,operatorLabel,site,spin);
		}

		assert(false);
		return 0;
	}

	void doBitmask()
	{
		SizeType n = geometry_.numberOfSites();
		bitmask_.resize(n);
		bitmask_[0]=1ul;
		for (SizeType i=1;i<n;i++)
			bitmask_[i] = bitmask_[i-1]<<1;
	}

	int doSign(WordType ket,SizeType i,SizeType j) const
	{
		throw PsimagLite::RuntimeError("BasisHeisenberg::doSign \n");
	}

	int getBraSzOrN(WordType& bra,
	                const WordType& ket1,
	                const WordType& ket2,
	                SizeType,
	                SizeType site,
	                SizeType spin) const
	{
		throw PsimagLite::RuntimeError("BasisHeisenberg::doSign \n");
	}

	SizeType logBase2(SizeType x) const
	{
		int ret = 0;
		while (x >>= 1) ++ret;
		return ret;
	}

	const GeometryType& geometry_;
	SizeType twiceS_;
	SizeType szPlusConst_;
	SizeType bits_;
	VectorWordType data_;
}; // class BasisHeisenberg

template<typename GeometryType>
std::ostream& operator<<(std::ostream& os,const BasisHeisenberg<GeometryType>& basis)
{
	for (SizeType i=0;i<basis.data_.size();i++)
		os<<i<<" "<<basis.data_[i]<<"\n";
	return os;
}

template<typename GeometryType>
typename BasisHeisenberg<GeometryType>::VectorWordType
BasisHeisenberg<GeometryType>::bitmask_;

} // namespace LanczosPlusPlus
#endif

