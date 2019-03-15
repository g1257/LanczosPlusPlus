/*
*/

#ifndef LANCZOS_BASIS_KITAEV_H
#define LANCZOS_BASIS_KITAEV_H

#include "BitManip.h"
#include "ProgramGlobals.h"
#include "../../Engine/BasisBase.h"

namespace LanczosPlusPlus {

template<typename GeometryType>
class BasisKitaev : public BasisBase<GeometryType> {

	typedef ProgramGlobals::PairIntType PairIntType;

public:

	static const SizeType TWICE_THE_SPIN = 1;
	static const SizeType BITS = 1;

	typedef BasisBase<GeometryType> BaseType;
	typedef typename BaseType::WordType WordType;
	typedef typename BaseType::VectorWordType VectorWordType;
	typedef typename BaseType::LabeledOperatorType LabeledOperatorType;

	static VectorWordType bitmask_;

	BasisKitaev(const GeometryType& geometry)
	    : geometry_(geometry), hilbert_(1)
	{
		SizeType sites = geometry_.numberOfSites();
		hilbert_ <<= sites;
		assert(bitmask_.size() == 0 || bitmask_.size() == sites);
		if (bitmask_.size()==0) doBitmask(sites);
	}

	PairIntType parts() const
	{
		throw PsimagLite::RuntimeError("BasisKitaev::parts()\n");
	}

	static const WordType& bitmask(SizeType i)
	{
		return bitmask_[i];
	}

	SizeType size() const { return hilbert_; }

	SizeType dofs() const
	{
		throw PsimagLite::RuntimeError("BasisKitaev::dofs() please check\n");
		return TWICE_THE_SPIN + 1;
	}

	SizeType hilbertOneSite(SizeType) const
	{
		return TWICE_THE_SPIN + 1;
	}

	SizeType perfectIndex(const VectorWordType&) const
	{
		throw PsimagLite::RuntimeError("BasisKitaev::perfectIndex kets\n");
	}

	SizeType perfectIndex(WordType ket,WordType) const
	{
		if (ket < hilbert_) return ket;

		throw PsimagLite::RuntimeError("perfectIndex: no index found\n");
	}

	SizeType perfectIndex(WordType,SizeType,SizeType) const
	{
		throw PsimagLite::RuntimeError("BasisKitaev::perfectIndex\n");
	}

	WordType operator()(SizeType i, SizeType) const
	{
		return i;
	}

	SizeType isThereAnElectronAt(WordType,
	                             WordType,
	                             SizeType,
	                             SizeType,
	                             SizeType) const
	{
		throw PsimagLite::RuntimeError("BasisKitaev::isThereAnElectronAt\n");
	}

	SizeType getN(WordType ket1,
	              WordType,
	              SizeType site,
	              SizeType,
	              SizeType) const
	{
		WordType mask = getMask();
		ket1 >>= (BITS*site);
		return (ket1 & mask);
	}

	int doSignGf(WordType, WordType,SizeType,SizeType,SizeType) const
	{
		throw PsimagLite::RuntimeError("BasisKitaev::doSignGf\n");
	}

	int doSign(WordType,
	           WordType,
	           SizeType,
	           SizeType,
	           SizeType,
	           SizeType,
	           SizeType) const
	{
		throw PsimagLite::RuntimeError("BasisKitaev::doSign\n");
	}

	PairIntType getBraIndex(WordType ket1,
	                        WordType ket2,
	                        const LabeledOperatorType&,
	                        SizeType site,
	                        SizeType spin,
	                        SizeType orb) const
	{
		throw PsimagLite::RuntimeError("BasisKitaev::getBraIndex() unimplemented yet\n");

		if (TWICE_THE_SPIN != 1)
			throw PsimagLite::RuntimeError("BasisKitaev::getBraIndex_ \n");

//		if (operatorLabel == ProgramGlobals::OPERATOR_SPLUS ||
//		        operatorLabel == ProgramGlobals::OPERATOR_SMINUS) {
//			return getBraIndexSplusSminus(ket1,ket2,operatorLabel,site,spin,orb);
//		}

//		return getBraIndex_(ket1,ket2,operatorLabel,site,spin,orb);
	}

	SizeType orbsPerSite(SizeType) const { return 1; }

	SizeType orbs() const { return 1; }

	void print(std::ostream& os, typename BaseType::PrintEnum binaryOrDecimal) const
	{
		SizeType hilbert = 1;
		hilbert <<= geometry_.numberOfSites();
		if (binaryOrDecimal == BaseType::PRINT_BINARY) {
			ProgramGlobals::printBasisBinary(os, hilbert, hilbert_);
		} else {
			ProgramGlobals::printBasisDecimal(os, 40, hilbert_);
		}
	}

	friend std::ostream& operator<<(std::ostream& os,
	                                const BasisKitaev& basis)
	{
		for (SizeType i=0;i<basis.data_.size();i++)
			os<<i<<" "<<basis.data_[i]<<"\n";
		return os;
	}

private:

	bool getBra(WordType& bra,
	            WordType ket,
	            WordType site1,
	            const LabeledOperatorType& val1,
	            SizeType site2,
	            SizeType val2) const
	{
		bra = ket;
		WordType mask1 = getMask();
		WordType mask2 = mask1;
		mask1 <<= (site1*BITS);
		bra &= (~mask1);

		mask2 <<= (site2*BITS);
		bra &= (~mask2);

		mask1 = val1.toUint();
		mask1 <<= (site1*BITS);
		bra |= mask1;

		mask2 = val2;
		mask2 <<= (site2*BITS);
		bra |= mask2;

		return true;
	}

	WordType getMask() const
	{
		SizeType mask = 1;
		for (SizeType i = 0; i < BITS; ++i)
			mask |= bitmask_[i];
		return mask;
	}

	static void doBitmask(SizeType n)
	{
		bitmask_.resize(n);
		bitmask_[0] = 1ul;
		for (SizeType i = 1; i < n; ++i)
			bitmask_[i] = bitmask_[i-1]<<1;
	}

	const GeometryType& geometry_;
	SizeType hilbert_;
}; // class BasisKitaev

template<typename GeometryType>
typename BasisKitaev<GeometryType>::VectorWordType
BasisKitaev<GeometryType>::bitmask_;

} // namespace LanczosPlusPlus
#endif

