/*
*/

#ifndef BASISRASHBASOC_H
#define BASISRASHBASOC_H
#include "../HubbardOneOrbital/BasisOneSpin.h"
#include "../../Engine/BasisBase.h"

namespace LanczosPlusPlus {

template<typename GeometryType>
class BasisRashbaSOC : public BasisBase<GeometryType> {

	enum {SPIN_UP = ProgramGlobals::SPIN_UP, SPIN_DOWN = ProgramGlobals::SPIN_DOWN};

public:

	typedef ProgramGlobals::PairIntType PairIntType;
	typedef BasisOneSpin BasisType;
	typedef BasisBase<GeometryType> BaseType;
	typedef typename BaseType::WordType WordType;
	typedef typename BaseType::VectorWordType VectorWordType;
	typedef typename BaseType::LabeledOperatorType LabeledOperatorType;
	static int const FERMION_SIGN = BasisType::FERMION_SIGN;

	BasisRashbaSOC(const GeometryType& geometry, SizeType ne)
	    : ne_(ne),
	      basis_(geometry.numberOfSites(), ne)
	{}

	PairIntType parts() const
	{
		throw PsimagLite::RuntimeError("parts() unimplemented\n");
	}

	static const WordType& bitmask(SizeType i)
	{
		return BasisType::bitmask(i);
	}

	SizeType size() const { return basis_.size(); }

	//! Spin up and spin down
	SizeType dofs() const
	{
		throw PsimagLite::RuntimeError("dofs() unimplemented\n");
	}

	SizeType hilbertOneSite(SizeType) const
	{
		return 4;
	}

	SizeType perfectIndex(const typename PsimagLite::Vector<WordType>::Type& kets) const
	{
		throw PsimagLite::RuntimeError("dofs() perfectIndex\n");
	}

	SizeType perfectIndex(WordType ket1,WordType ket2) const
	{
		throw PsimagLite::RuntimeError("dofs() perfectIndex\n");
	}

	SizeType perfectIndex(WordType newKet,
	                      SizeType ispace,
	                      SizeType spinOfNew) const
	{
		throw PsimagLite::RuntimeError("dofs() perfectIndex\n");
	}

	WordType operator()(SizeType i,SizeType spin) const
	{
		throw PsimagLite::RuntimeError("operator()\n");
	}

	PairIntType getBraIndex(WordType ket1,
	                                WordType ket2,
	                                const LabeledOperatorType&,
	                                SizeType site,
	                                SizeType spin,
	                                SizeType orb) const
	{
		throw PsimagLite::RuntimeError("getBraIndex()\n");
	}

	int doSign(WordType ket1,
	                   WordType ket2,
	                   SizeType i,
	                   SizeType orb1,
	                   SizeType j,
	                   SizeType orb2,
	                   SizeType spin) const
	{
		throw PsimagLite::RuntimeError("doSign()\n");
	}

	int doSignGf(WordType a,
	                     WordType b,
	                     SizeType ind,
	                     SizeType spin,
	                     SizeType orb) const
	{
		throw PsimagLite::RuntimeError("doSignGf()\n");
	}

	int doSignSpSm(WordType a,
                             WordType b,
                             SizeType ind,
                             SizeType spin,
                             SizeType orb) const
	{
		return 1;
	}

	SizeType isThereAnElectronAt(WordType ket1,
	                                     WordType ket2,
	                                     SizeType site,
	                                     SizeType spin,
	                                     SizeType orb) const
	{
		throw PsimagLite::RuntimeError("isThereAnElectronAt()\n");
	}

	SizeType orbsPerSite(SizeType i) const  { return 1; }

	SizeType orbs() const { return 1; }

	SizeType getN(WordType ket1,
	                      WordType ket2,
	                      SizeType site,
	                      SizeType spin,
	                      SizeType orb) const
	{
		throw PsimagLite::RuntimeError("getN()\n");
	}

	bool getBra(WordType&,
	                    WordType,
	                    WordType,
	                    const LabeledOperatorType&,
	                    SizeType,
	                    SizeType) const
	{
		throw PsimagLite::RuntimeError("getBra()\n");
	}

	void print(std::ostream&, typename BaseType::PrintEnum) const
	{
		throw PsimagLite::RuntimeError("print()\n");
	}

private:

	SizeType ne_;
	BasisType basis_;

}; // class BasisRashbaSOC

} // namespace LanczosPlusPlus
#endif

