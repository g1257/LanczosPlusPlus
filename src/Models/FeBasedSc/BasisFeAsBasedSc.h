
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

#ifndef BASIS_FEASBASED_SC_H
#define BASIS_FEASBASED_SC_H
#include "BasisOneSpinFeAs.h"
#include "../../Engine/BasisBase.h"

namespace LanczosPlusPlus {

template<typename GeometryType>
class BasisFeAsBasedSc : public BasisBase<GeometryType> {

	typedef ProgramGlobals::PairIntType PairIntType;

	static SizeType orbitals_;
	static const SizeType OPERATOR_C = ProgramGlobals::OPERATOR_C;
	static const SizeType OPERATOR_CDAGGER = ProgramGlobals::OPERATOR_CDAGGER;
	static const SizeType SPIN_UP = ProgramGlobals::SPIN_UP;

public:

	typedef BasisBase<GeometryType> BaseType;
	typedef typename BaseType::WordType WordType;
	typedef typename BaseType::VectorWordType VectorWordType;
	typedef BasisOneSpinFeAs BasisType;
	static int const FERMION_SIGN = BasisType::FERMION_SIGN;

	BasisFeAsBasedSc(const GeometryType& geometry, SizeType nup,SizeType ndown,SizeType orbitals)
	: basis1_(geometry.numberOfSites(),nup,orbitals),
	  basis2_(geometry.numberOfSites(),ndown,orbitals)
	{
		orbitals_ = orbitals;
	}

	BasisFeAsBasedSc(const GeometryType& geometry, SizeType nup,SizeType ndown)
	: basis1_(geometry.numberOfSites(),nup,orbitals_),
	  basis2_(geometry.numberOfSites(),ndown,orbitals_)
	{}

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
		return basis1_.size()*basis2_.size();
	}

	WordType operator()(SizeType i,SizeType spin) const
	{
		SizeType y = i/basis1_.size();
		SizeType x = i%basis1_.size();
		return (spin==SPIN_UP) ? basis1_[x] : basis2_[y];
	}

	SizeType perfectIndex(const VectorWordType& kets) const
	{
		assert(kets.size()==2);
		return  perfectIndex(kets[0],kets[1]);
	}

	SizeType perfectIndex(WordType ket1,WordType ket2) const
	{
		return basis1_.perfectIndex(ket1) + basis2_.perfectIndex(ket2)*basis1_.size();
	}

	SizeType perfectIndex(WordType newKet,
	                      SizeType ispace,
	                      SizeType spinOfNew) const
	{
		throw PsimagLite::RuntimeError("perfectIndex\n");
	}

	SizeType getN(SizeType i,SizeType spin,SizeType orb) const
	{
		SizeType y = i/basis1_.size();
		SizeType x = i%basis1_.size();
		return (spin==SPIN_UP) ? basis1_.getN(x,orb) : basis2_.getN(y,orb);
	}

	SizeType getN(WordType ket1, WordType ket2, SizeType site,SizeType spin,SizeType orb) const
	{
		return (spin==SPIN_UP) ? basis1_.getN(ket1,site,orb) : basis2_.getN(ket2,site,orb);
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
			return basis1_.doSign(ket1,i,orb1,j,orb2);
		}
		return basis2_.doSign(ket2,i,orb1,j,orb2);
	}

	int doSignGf(WordType a, WordType b,SizeType ind,SizeType spin,SizeType orb) const
	{
		if (spin==SPIN_UP) return basis1_.doSignGf(a,ind,orb);

		int s=(PsimagLite::BitManip::count(a) & 1) ? -1 : 1; // Parity of up
		int s2 = basis2_.doSignGf(b,ind,orb);

		return s*s2;
	}

	SizeType isThereAnElectronAt(WordType ket1,
	                             WordType ket2,
	                             SizeType site,
	                             SizeType spin,
	                             SizeType orb) const
	{
		if (spin==SPIN_UP)
			return basis1_.isThereAnElectronAt(ket1,site,orb);
		return basis2_.isThereAnElectronAt(ket2,site,orb);
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

	SizeType orbsPerSite(SizeType i) const
	{
		throw PsimagLite::RuntimeError("orbsPerSite\n");
	}

	SizeType orbs() const
	{
		throw PsimagLite::RuntimeError("orbs\n");
	}

private:

	int getBraIndexCorCdaggerOrN(WordType ket1,
	                             WordType ket2,
	                             SizeType what,
	                             SizeType site,
	                             SizeType spin,
	                             SizeType orb) const
	{

		WordType bra  =0;
		bool b = getBraCorCdaggerOrN(bra,ket1,ket2,what,site,spin,orb);
		if (!b) return -1;
		return (spin==SPIN_UP) ? perfectIndex(bra,ket2) :
		       perfectIndex(ket1,bra);
	}

	int getBraIndexSplusOrSminus(WordType ket1,
	                             WordType ket2,
	                             SizeType what,
	                             SizeType site,
	                             SizeType orb) const
	{

		WordType bra1  =0;
		WordType bra2  =0;
		bool b = getBraSplusOrSminus(bra1,bra2,ket1,ket2,what,site,orb);
		if (!b) return -1;
		return perfectIndex(bra1,bra2);
	}

	bool hasNewPartsCorCdagger(std::pair<SizeType,SizeType>& newParts,
	                           SizeType what,
	                           SizeType spin,
	                           const std::pair<SizeType,SizeType>& orbs) const
	{
		int newPart1=basis1_.electrons();
		int newPart2=basis2_.electrons();

		if (spin==SPIN_UP) newPart1 = basis1_.newPartCorCdagger(what,orbs.first);
		else newPart2 = basis2_.newPartCorCdagger(what,orbs.second);

		if (newPart1<0 || newPart2<0) return false;

		if (newPart1==0 && newPart2==0) return false;
		newParts.first = SizeType(newPart1);
		newParts.second = SizeType(newPart2);
		return true;
	}

	bool hasNewPartsSplusOrSminus(std::pair<SizeType,SizeType>& newParts,
	                              SizeType what,
	                              const std::pair<SizeType,SizeType>& orbs) const
	{
		int c1 = (what==ProgramGlobals::OPERATOR_SPLUS) ? 1 : -1;
		int c2 = (what==ProgramGlobals::OPERATOR_SPLUS) ? -1 : 1;

		int newPart1 = basis1_.hasNewPartsSplusOrSminus(c1,orbs.first);
		int newPart2 = basis2_.hasNewPartsSplusOrSminus(c2,orbs.second);

		if (newPart1<0 || newPart2<0) return false;

		if (newPart1==0 && newPart2==0) return false;
		newParts.first = SizeType(newPart1);
		newParts.second = SizeType(newPart2);
		return true;
	}

	bool getBraCorCdaggerOrN(WordType& bra,
	                         const WordType& ket1,
	                         const WordType& ket2,
	                         SizeType what,
	                         SizeType site,
	                         SizeType spin,
	                         SizeType orb) const
	{
		return (spin==SPIN_UP) ? basis1_.getBra(bra,ket1,what,site,orb) :
		       basis2_.getBra(bra,ket2,what,site,orb);
	}

	bool getBraSplusOrSminus(WordType& bra1,
	                         WordType& bra2,
	                         const WordType& ket1,
	                         const WordType& ket2,
	                         SizeType what,
	                         SizeType site,
	                         SizeType orb) const
	{
		SizeType what1 = (what==ProgramGlobals::OPERATOR_SPLUS) ? OPERATOR_CDAGGER : OPERATOR_C;
		SizeType what2 = (what==ProgramGlobals::OPERATOR_SPLUS) ? OPERATOR_C : OPERATOR_CDAGGER;
		bool b1 = basis1_.getBra(bra1,ket1,what1,site,orb);
		bool b2 = basis2_.getBra(bra2,ket2,what2,site,orb);
		return (b1 & b2);
	}

	BasisType basis1_,basis2_;
}; // class BasisFeAsBasedSc

template<typename GeometryType>
SizeType BasisFeAsBasedSc<GeometryType>::orbitals_=2;

} // namespace LanczosPlusPlus
#endif

