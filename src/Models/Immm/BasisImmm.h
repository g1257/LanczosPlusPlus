
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

#ifndef BASIS_IMMM_H
#define BASIS_IMMM_H
#include "BasisOneSpinImmm.h"
#include "ProgramGlobals.h"
#include "BasisBase.h"

namespace LanczosPlusPlus {

template<typename GeometryType>
class BasisImmm : public BasisBase<GeometryType> {

	typedef ProgramGlobals::PairIntType PairIntType;

public:

	typedef BasisBase<GeometryType> BaseType;
        typedef typename BaseType::WordType WordType;
        typedef typename BaseType::VectorWordType VectorWordType;
	typedef BasisOneSpinImmm BasisType;

	class OrbsPerSite : public PsimagLite::Vector<SizeType>::Type {

		public:

			OrbsPerSite(const GeometryType& geometry)
			: PsimagLite::Vector<SizeType>::Type(geometry.numberOfSites())
			{
				typename GeometryType::AdditionalDataType additional;
				additional.type1 = 0;
				additional.TYPE_C = 0;
				for (SizeType i=0;i<this->size();i++) {
					geometry.fillAdditionalData(additional,0,i,0);
					this->operator[](i) = (additional.type1==additional.TYPE_C) ? 1 : 2;
				}
			}
		}; // class OrbsPerSite

		static int const FERMION_SIGN = BasisType::FERMION_SIGN;

		BasisImmm(const GeometryType& geometry, SizeType nup,SizeType ndown)
		: orbsPerSite_(geometry),
		  basis1_(orbsPerSite_,nup),
		  basis2_(orbsPerSite_,ndown)
		{}

		static const WordType& bitmask(SizeType i)
		{
			return BasisType::bitmask(i);
		}

		SizeType dofs() const
		{
			throw std::runtime_error("Wrong way!\n");
		}

		SizeType size() const { return basis1_.size()*basis2_.size(); }

		WordType operator()(SizeType i,SizeType spin) const
		{
			SizeType y = i/basis1_.size();
			SizeType x = i%basis1_.size();
			return (spin==ProgramGlobals::SPIN_UP) ? basis1_[x] : basis2_[y];
		}

		SizeType perfectIndex(const VectorWordType& ket1) const
		{
			throw std::runtime_error("Wrong way!\n");
		}

		SizeType perfectIndex(WordType newKet,SizeType ispace,SizeType spinOfNew) const
		{
			if (spinOfNew==ProgramGlobals::SPIN_UP) {
				SizeType oldIndex1 = ispace / basis1_.size();
				return basis1_.perfectIndex(newKet) + oldIndex1*basis1_.size();
			}
			SizeType oldIndex2 = ispace % basis2_.size();
			return oldIndex2 + basis2_.perfectIndex(newKet) * basis2_.size();
		}

		SizeType perfectIndex(WordType ket1,WordType ket2) const
		{
			return basis1_.perfectIndex(ket1) + basis2_.perfectIndex(ket2)*basis1_.size();
		}

		SizeType getN(SizeType i,SizeType spin,SizeType orb) const
		{
			SizeType y = i/basis1_.size();
			SizeType x = i%basis1_.size();
			return (spin==ProgramGlobals::SPIN_UP) ? basis1_.getN(x,orb) : basis2_.getN(y,orb);
		}

		SizeType getN(WordType ket,SizeType site,SizeType spin,SizeType orb) const
		{
			return (spin==ProgramGlobals::SPIN_UP) ? basis1_.getN(ket,site,orb) : basis2_.getN(ket,site,orb);
		}

		PairIntType getBraIndex(WordType ket1,
		                        WordType ket2,
		                        SizeType what,
		                        SizeType site,
		                        SizeType spin,
		                        SizeType orb) const
		{
			WordType bra = 0;
			bool b = getBra(bra,ket1,ket2,what,site,spin,orb);
			if (!b) return PairIntType(-1,1);
			int tmp = (spin==ProgramGlobals::SPIN_UP) ? perfectIndex(bra,ket2) : perfectIndex(ket1,bra);
			return PairIntType(tmp,1);
		}

		int doSign(SizeType i,SizeType site,SizeType sector) const
		{
			SizeType y = i/basis1_.size();
			SizeType x = i%basis1_.size();
			SizeType spin = sector/2;
			SizeType orb = (sector & 1);
			if (spin==ProgramGlobals::SPIN_UP) {
				return basis1_.doSign(x,site,orb);
			}
			SizeType c = basis1_.getN(x);
			int ret = 1;
			if (c&1) ret = FERMION_SIGN;
			return ret * basis2_.doSign(y,site,orb);
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
			if (spin==ProgramGlobals::SPIN_UP) {
				return basis1_.doSign(ket1,i,orb1,j,orb2);
			}
			return basis2_.doSign(ket2,i,orb1,j,orb2);
		}

		int doSignGf(WordType a, WordType b,SizeType ind,SizeType spin,SizeType orb) const
		{
			if (spin==ProgramGlobals::SPIN_UP) {
				return basis1_.doSignGf(a,ind,orb);
			}
			int s=(PsimagLite::BitManip::count(a) & 1) ? -1 : 1; // Parity of up
			return s*basis2_.doSignGf(b,ind,orb);
		}

		SizeType isThereAnElectronAt(SizeType ket1,
		                           SizeType ket2,
		                           SizeType site,
		                           SizeType spin,
		                           SizeType orb) const
		{
			if (spin==ProgramGlobals::SPIN_UP)
				return basis1_.isThereAnElectronAt(ket1,site,orb);
			return basis2_.isThereAnElectronAt(ket2,site,orb);
		}

		SizeType orbsPerSite(SizeType i) const { return orbsPerSite_[i]; }

		SizeType orbs() const { return orbsPerSite_[0]; }

		SizeType electrons(SizeType what) const
		{
			return (what==ProgramGlobals::SPIN_UP) ? basis1_.electrons() : basis2_.electrons();
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

	private:

		bool hasNewPartsCorCdagger(std::pair<SizeType,SizeType>& newParts,
		                           SizeType what,
		                           SizeType spin,
		                           const std::pair<SizeType,SizeType>& orbs) const
		{
			int newPart1=basis1_.electrons();
			int newPart2=basis2_.electrons();

			if (spin==ProgramGlobals::SPIN_UP) newPart1 = basis1_.newPartCorCdagger(what,orbs.first);
			else newPart2 = basis2_.newPartCorCdagger(what,orbs.second);

			if (newPart1<0 || newPart2<0) return false;

			if (newPart1==0 && newPart2==0) return false;
			newParts.first = SizeType(newPart1);
			newParts.second = SizeType(newPart2);
			return true;
		}

		bool getBra(WordType& bra,
					const WordType& ket1,
					const WordType& ket2,
					SizeType what,
					SizeType site,
					SizeType spin,
					SizeType orb) const
		{
			return (spin==ProgramGlobals::SPIN_UP) ? basis1_.getBra(bra,ket1,what,site,orb) :
									 basis2_.getBra(bra,ket2,what,site,orb);
		}

		OrbsPerSite orbsPerSite_;
		BasisType basis1_,basis2_;

	}; // class BasisImmm
} // namespace LanczosPlusPlus
#endif
