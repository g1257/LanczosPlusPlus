/*
Copyright (c) 2009-2016, UT-Battelle, LLC
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

#ifndef BASIS_FEASBASED_SPINORBIT_H
#define BASIS_FEASBASED_SPINORBIT_H
#include "../../Engine/BasisBase.h"

namespace LanczosPlusPlus {

template<typename GeometryType_>
class BasisFeAsSpinOrbit : public BasisBase<GeometryType_> {

	typedef ProgramGlobals::PairIntType PairIntType;

	static SizeType orbitals_;
	static const SizeType OPERATOR_C = ProgramGlobals::OPERATOR_C;
	static const SizeType OPERATOR_CDAGGER = ProgramGlobals::OPERATOR_CDAGGER;
	static const SizeType SPIN_UP = ProgramGlobals::SPIN_UP;

public:

	typedef GeometryType_ GeometryType;
	typedef BasisBase<GeometryType> BaseType;
	typedef typename BaseType::WordType WordType;
	typedef typename BaseType::VectorWordType VectorWordType;
	typedef BasisOneSpinFeAs BasisType;
	static int const FERMION_SIGN = BasisType::FERMION_SIGN;

	BasisFeAsSpinOrbit(const GeometryType& geometry,
	                   SizeType nup,
	                   SizeType ndown,
	                   SizeType orbitals)
	{
		orbitals_ = orbitals;
	}

	static const WordType& bitmask(SizeType i)
	{
		return BasisType::bitmask(i);
	}

	SizeType dofs() const
	{
		return 2*orbitals_;
	}

	SizeType size() const
	{
		return 0;
	}

	virtual SizeType hilbertOneSite(SizeType) const
	{
		return (1<<(orbitals_*2));
	}

	WordType operator()(SizeType i,SizeType spin) const
	{
		return (spin==SPIN_UP) ? 0 : 0;
	}

	SizeType perfectIndex(const VectorWordType& kets) const
	{
		assert(kets.size()==2);
		return  perfectIndex(kets[0],kets[1]);
	}

	SizeType perfectIndex(WordType ket1,WordType ket2) const
	{
		return 0;
	}

	SizeType perfectIndex(WordType,
	                      SizeType,
	                      SizeType) const
	{
		throw PsimagLite::RuntimeError("perfectIndex\n");
	}

	SizeType getN(SizeType i,SizeType spin,SizeType orb) const
	{
		return (spin==SPIN_UP) ? 0 : 0;
	}

	SizeType getN(WordType ket1, WordType ket2, SizeType site,SizeType spin,SizeType orb) const
	{
		return (spin==SPIN_UP) ? 0 : 0;
	}

	PairIntType getBraIndex(WordType ket1,
	                        WordType ket2,
	                        SizeType what,
	                        SizeType site,
	                        SizeType spin,
	                        SizeType orb) const
	{
		if (what==ProgramGlobals::OPERATOR_C ||
		    what==ProgramGlobals::OPERATOR_CDAGGER ||
		    what==ProgramGlobals::OPERATOR_N)
		  return PairIntType(getBraIndexCorCdaggerOrN(ket1,ket2,what,site,spin,orb),1);
		if (what==ProgramGlobals::OPERATOR_SPLUS || what==ProgramGlobals::OPERATOR_SMINUS)
			return PairIntType(getBraIndexSplusOrSminus(ket1,ket2,what,site,orb),1);
		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) +  "\n";
		str += PsimagLite::String("getBraIndex: unsupported operator ");
		str += ProgramGlobals::id2Operator(what) + "\n";
		throw std::runtime_error(str.c_str());
	}

	int doSign(WordType ket1,
	           WordType ket2,
	           SizeType i,
	           SizeType orb1,
	           SizeType j,
	           SizeType orb2,
	           SizeType spin) const
	{
		if (i > j) {
			std::cerr<<"FATAL: At doSign\n";
			std::cerr<<"INFO: i="<<i<<" j="<<j<<std::endl;
			std::cerr<<"AT: "<<__FILE__<<" : "<<__LINE__<<std::endl;
			throw std::runtime_error("FeBasedSc::doSign(...)\n");
		}

		if (spin==SPIN_UP) {
			return 0;
		}
		return 0;
	}

	int doSignGf(WordType a, WordType b,SizeType ind,SizeType spin,SizeType orb) const
	{
		if (spin==SPIN_UP) return 0;

		int s=(PsimagLite::BitManip::count(a) & 1) ? -1 : 1; // Parity of up
		int s2 = 0;

		return s*s2;
	}

	int doSignSpinOrbit(WordType a,
	                    WordType b,
	                    SizeType ind,
	                    SizeType spin1,
	                    SizeType orb1,
	                    SizeType spin2,
	                    SizeType orb2) const
	{
		return 0;
	}

	int doSignSpSm(WordType a, WordType b,SizeType ind,SizeType spin,SizeType orb) const
	{
		return 0;
	}

	SizeType isThereAnElectronAt(WordType ket1,
	                             WordType ket2,
	                             SizeType site,
	                             SizeType spin,
	                             SizeType orb) const
	{
		return 0;
	}

	bool hasNewParts(std::pair<SizeType,SizeType>& newParts,
	                 SizeType what,
	                 SizeType spin,
	                 const std::pair<SizeType,SizeType>& orbs) const
	{
		if (what==ProgramGlobals::OPERATOR_C || what==ProgramGlobals::OPERATOR_CDAGGER)
			return hasNewPartsCorCdagger(newParts,what,spin,orbs);
		if (what==ProgramGlobals::OPERATOR_SPLUS || what==ProgramGlobals::OPERATOR_SMINUS)
			return hasNewPartsSplusOrSminus(newParts,what,orbs);
		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) +  "\n";
		str += PsimagLite::String("hasNewParts: unsupported operator ");
		str += ProgramGlobals::id2Operator(what) + "\n";
		throw std::runtime_error(str.c_str());
	}

	SizeType orbsPerSite(SizeType) const
	{
		throw PsimagLite::RuntimeError("orbsPerSite\n");
	}

	SizeType orbs() const
	{
		throw PsimagLite::RuntimeError("orbs\n");
	}

	void print(std::ostream& os, typename BaseType::PrintEnum binaryOrDecimal) const
	{
//		bool isBinary = (binaryOrDecimal == BaseType::PRINT_BINARY);
//		os<<"\tUp sector\n";
//		basis1_.print(os,isBinary);
//		os<<"\tDown sector\n";
//		basis2_.print(os,isBinary);
	}

	virtual bool getBra(WordType&,
	                    WordType,
	                    WordType,
	                    SizeType,
	                    SizeType,
	                    SizeType) const
	{
		throw PsimagLite::RuntimeError("BasisFeAs: getBra\n");
	}

private:

	int getBraIndexCorCdaggerOrN(WordType ket1,
	                             WordType ket2,
	                             SizeType what,
	                             SizeType site,
	                             SizeType spin,
	                             SizeType orb) const
	{
		return 0;
	}

	int getBraIndexSplusOrSminus(WordType ket1,
	                             WordType ket2,
	                             SizeType what,
	                             SizeType site,
	                             SizeType orb) const
	{

		return 0;
	}

	bool hasNewPartsCorCdagger(std::pair<SizeType,SizeType>& newParts,
	                           SizeType what,
	                           SizeType spin,
	                           const std::pair<SizeType,SizeType>& orbs) const
	{
		return false;
	}

	bool hasNewPartsSplusOrSminus(std::pair<SizeType,SizeType>& newParts,
	                              SizeType what,
	                              const std::pair<SizeType,SizeType>& orbs) const
	{
		return false;
	}

	bool getBraCorCdaggerOrN(WordType& bra,
	                         const WordType& ket1,
	                         const WordType& ket2,
	                         SizeType what,
	                         SizeType site,
	                         SizeType spin,
	                         SizeType orb) const
	{
		return false;
	}

	bool getBraSplusOrSminus(WordType& bra1,
	                         WordType& bra2,
	                         const WordType& ket1,
	                         const WordType& ket2,
	                         SizeType what,
	                         SizeType site,
	                         SizeType orb) const
	{
		return false;
	}
}; // class BasisFeAsSpinOrbit

template<typename GeometryType>
SizeType BasisFeAsSpinOrbit<GeometryType>::orbitals_=2;

} // namespace LanczosPlusPlus
#endif

