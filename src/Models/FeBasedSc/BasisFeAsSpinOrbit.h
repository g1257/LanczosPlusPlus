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
#include "BasisOneSpinFeAs.h"

namespace LanczosPlusPlus {

template<typename GeometryType_>
class BasisFeAsSpinOrbit : public BasisBase<GeometryType_> {

	typedef ProgramGlobals::PairIntType PairIntType;

	static SizeType orbitals_;
	static const SizeType OPERATOR_C = ProgramGlobals::OPERATOR_C;
	static const SizeType OPERATOR_CDAGGER = ProgramGlobals::OPERATOR_CDAGGER;
	static const SizeType SPIN_UP = ProgramGlobals::SPIN_UP;

public:

	typedef BasisOneSpinFeAs BasisOneSpinType;
	typedef GeometryType_ GeometryType;
	typedef BasisBase<GeometryType> BaseType;
	typedef typename BaseType::WordType WordType;
	typedef std::pair<WordType,WordType> PairWordType;
	typedef typename BaseType::VectorWordType VectorWordType;
	typedef typename PsimagLite::Vector<PairWordType>::Type VectorPairWordType;
	static int const FERMION_SIGN = -1;

	BasisFeAsSpinOrbit(const GeometryType& geometry,
	                   SizeType nup1,
	                   SizeType ndown1,
	                   SizeType orbitals)
	    : geometry_(geometry),
	      basis1_(geometry.numberOfSites(),1,orbitals) // bogus, just to use some functions
	{
		orbitals_ = orbitals;
		SizeType ne = nup1 + ndown1;
		for (SizeType nup = 0; nup <= ne; ++nup) {
			SizeType ndown = ne - nup;
			BasisOneSpinType basis1(geometry.numberOfSites(),nup,orbitals);
			BasisOneSpinType basis2(geometry.numberOfSites(),ndown,orbitals);
			for (SizeType i = 0; i < basis1.size(); ++i)
				for (SizeType j = 0; j < basis2.size(); ++j)
					basis_.push_back(PairWordType(basis1[i],basis2[j]));
		}

	}

	static const WordType& bitmask(SizeType i)
	{
		return BasisOneSpinType::bitmask(i);
	}

	SizeType dofs() const
	{
		return 2*orbitals_;
	}

	SizeType size() const
	{
		return basis_.size();
	}

	virtual SizeType hilbertOneSite(SizeType) const
	{
		return (1<<(orbitals_*2));
	}

	WordType operator()(SizeType i,SizeType spin) const
	{
		return (spin==SPIN_UP) ? basis_[i].first : basis_[i].second;
	}

	SizeType perfectIndex(const VectorWordType& kets) const
	{
		assert(kets.size()==2);
		return  perfectIndex(kets[0],kets[1]);
	}

	SizeType perfectIndex(WordType ket1,WordType ket2) const
	{
		PairWordType combined(ket1,ket2);
		typename VectorPairWordType::const_iterator it = std::find(basis_.begin(),
		                                                           basis_.end(),
		                                                           combined);
		assert(it != basis_.end());
		return it - basis_.begin();
	}

	SizeType perfectIndex(WordType,
	                      SizeType,
	                      SizeType) const
	{
		throw PsimagLite::RuntimeError("perfectIndex\n");
	}

	SizeType getN(WordType ket1,
	              WordType ket2,
	              SizeType site,
	              SizeType spin,
	              SizeType orb) const
	{
		return (spin==SPIN_UP) ? basis1_.getN(ket1,site,orb) :
		                         basis1_.getN(ket2,site,orb);
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

		return (spin==SPIN_UP) ? basis1_.doSign(ket1,i,orb1,j,orb2) :
		                         basis1_.doSign(ket2,i,orb1,j,orb2);
	}

	int doSignGf(WordType a,
	             WordType b,
	             SizeType ind,
	             SizeType spin,
	             SizeType orb) const
	{
		if (spin==SPIN_UP) return basis1_.doSignGf(a,ind,orb);

		// Parity of up
		int s=(PsimagLite::BitManip::count(a) & 1) ? -1 : 1;
		int s2 = basis1_.doSignGf(b,ind,orb);

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
		if (spin1==spin2) {
			return (spin1 == SPIN_UP) ? basis1_.doSign(a,ind,orb1,orb2) :
			                            basis1_.doSign(b,ind,orb1,orb2);
		}

		int x = (spin1 == SPIN_UP) ? 1 : -1;
		// Parity of up
		int s=(PsimagLite::BitManip::count(a) & 1) ? -1 : 1;
		if (spin1 == SPIN_UP) return x*s*basis1_.doSignGf(a,ind,orb1)*
		        basis1_.doSignGf(b,ind,orb2);

		return x*s*basis1_.doSignGf(a,ind,orb2)*basis1_.doSignGf(b,ind,orb1);
	}

	int doSignSpSm(WordType a,
	               WordType b,
	               SizeType ind,
	               SizeType spin,
	               SizeType orb) const
	{
		if (spin==SPIN_UP) { // spin here means S^\dagger
			// FIXME: Count over a (up)
			return basis1_.doSignGf(a,ind,orb)*basis1_.doSignGf(b,ind,orb);
		}

		// FIXME: Count over a + 1
		return basis1_.doSignGf(a,ind,orb)*basis1_.doSignGf(b,ind,orb);
	}

	SizeType isThereAnElectronAt(WordType ket1,
	                             WordType ket2,
	                             SizeType site,
	                             SizeType spin,
	                             SizeType orb) const
	{
		return (spin == SPIN_UP) ? basis1_.isThereAnElectronAt(ket1,site,orb) :
		                           basis1_.isThereAnElectronAt(ket2,site,orb);
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

	void print(std::ostream& os,
	           typename BaseType::PrintEnum binaryOrDecimal) const
	{
		bool isBinary = (binaryOrDecimal == BaseType::PRINT_BINARY);
		SizeType hilbert = 1;
		hilbert <<= (orbitals_*geometry_.numberOfSites());
		if (isBinary) {
			ProgramGlobals::printBasisBinary(os,hilbert,basis_);
		} else {
			ProgramGlobals::printBasisDecimal(os,40,basis_);
		}
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
		throw PsimagLite::RuntimeError("UNIMPLEMENTED: hasNewPartsCorCdagger\n");
	}

	bool hasNewPartsSplusOrSminus(std::pair<SizeType,SizeType>& newParts,
	                              SizeType what,
	                              const std::pair<SizeType,SizeType>& orbs) const
	{
		throw PsimagLite::RuntimeError("UNIMPLEMENTED: hasNewPartsSplusOrSminus\n");
	}

	bool getBraCorCdaggerOrN(WordType& bra,
	                         const WordType& ket1,
	                         const WordType& ket2,
	                         SizeType what,
	                         SizeType site,
	                         SizeType spin,
	                         SizeType orb) const
	{
		throw PsimagLite::RuntimeError("UNIMPLEMENTED: getBraCorCdaggerOrN\n");
	}

	bool getBraSplusOrSminus(WordType& bra1,
	                         WordType& bra2,
	                         const WordType& ket1,
	                         const WordType& ket2,
	                         SizeType what,
	                         SizeType site,
	                         SizeType orb) const
	{
		throw PsimagLite::RuntimeError("UNIMPLEMENTED: getBraSplusOrSminus\n");
	}


	const GeometryType& geometry_;
	VectorPairWordType basis_;
	BasisOneSpinType basis1_;
}; // class BasisFeAsSpinOrbit

template<typename GeometryType>
SizeType BasisFeAsSpinOrbit<GeometryType>::orbitals_=2;

} // namespace LanczosPlusPlus
#endif

