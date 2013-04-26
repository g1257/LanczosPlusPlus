
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

namespace LanczosPlusPlus {

	template<typename GeometryType>
	class BasisImmm {
	public:

		typedef BasisOneSpinImmm BasisType;
		typedef BasisType::WordType WordType;

		class OrbsPerSite : public std::vector<size_t> {

		public:

			OrbsPerSite(const GeometryType& geometry)
			: std::vector<size_t>(geometry.numberOfSites())
			{
				typename GeometryType::AdditionalDataType additional;
				additional.type1 = 0;
				additional.TYPE_C = 0;
				for (size_t i=0;i<this->size();i++) {
					geometry.fillAdditionalData(additional,0,i,0);
					this->operator[](i) = (additional.type1==additional.TYPE_C) ? 1 : 2;
				}
			}
		}; // class OrbsPerSite

		static int const FERMION_SIGN = BasisType::FERMION_SIGN;

		BasisImmm(const GeometryType& geometry, size_t nup,size_t ndown)
		: orbsPerSite_(geometry),
		  basis1_(orbsPerSite_,nup),
		  basis2_(orbsPerSite_,ndown)
		{}

		static const WordType& bitmask(size_t i)
		{
			return BasisType::bitmask(i);
		}

		size_t dofs() const
		{
			throw std::runtime_error("Wrong way!\n");
		}

		size_t size() const { return basis1_.size()*basis2_.size(); }

		const WordType& operator()(size_t i,size_t spin) const
		{
			size_t y = i/basis1_.size();
			size_t x = i%basis1_.size();
			return (spin==ProgramGlobals::SPIN_UP) ? basis1_[x] : basis2_[y];
		}

		size_t perfectIndex(const std::vector<WordType>& ket1) const
		{
			throw std::runtime_error("Wrong way!\n");
		}

		size_t perfectIndex(WordType newKet,size_t ispace,size_t spinOfNew) const
		{
			if (spinOfNew==ProgramGlobals::SPIN_UP) {
				size_t oldIndex1 = ispace / basis1_.size();
				return basis1_.perfectIndex(newKet) + oldIndex1*basis1_.size();
			}
			size_t oldIndex2 = ispace % basis2_.size();
			return oldIndex2 + basis2_.perfectIndex(newKet) * basis2_.size();
		}

		size_t perfectIndex(WordType ket1,WordType ket2) const
		{
			return basis1_.perfectIndex(ket1) + basis2_.perfectIndex(ket2)*basis1_.size();
		}

		size_t getN(size_t i,size_t spin,size_t orb) const
		{
			size_t y = i/basis1_.size();
			size_t x = i%basis1_.size();
			return (spin==ProgramGlobals::SPIN_UP) ? basis1_.getN(x,orb) : basis2_.getN(y,orb);
		}

		size_t getN(WordType ket,size_t site,size_t spin,size_t orb) const
		{
			return (spin==ProgramGlobals::SPIN_UP) ? basis1_.getN(ket,site,orb) : basis2_.getN(ket,site,orb);
		}

		size_t getBraIndex(const WordType& ket1,
		                   const WordType& ket2,
		                   size_t what,
		                   size_t site,
		                   size_t spin,
		                   size_t orb) const
		{
			WordType bra = 0;
			bool b = getBra(bra,ket1,ket2,what,site,spin,orb);
			if (!b) return -1;
			return (spin==ProgramGlobals::SPIN_UP) ? perfectIndex(bra,ket2) : perfectIndex(ket1,bra);
		}

		int doSign(size_t i,size_t site,size_t sector) const
		{
			size_t y = i/basis1_.size();
			size_t x = i%basis1_.size();
			size_t spin = sector/2;
			size_t orb = (sector & 1);
			if (spin==ProgramGlobals::SPIN_UP) {
				return basis1_.doSign(x,site,orb);
			}
			size_t c = basis1_.getN(x);
			int ret = 1;
			if (c&1) ret = FERMION_SIGN;
			return ret * basis2_.doSign(y,site,orb);
		}

		int doSign(WordType ket1,
		           WordType ket2,
		           size_t i,
		           size_t orb1,
		           size_t j,
		           size_t orb2,
		           size_t spin) const
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

		int doSignGf(WordType a, WordType b,size_t ind,size_t spin,size_t orb) const
		{
			if (spin==ProgramGlobals::SPIN_UP) {
				return basis1_.doSignGf(a,ind,orb);
			}
			int s=(PsimagLite::BitManip::count(a) & 1) ? -1 : 1; // Parity of up
			return s*basis2_.doSignGf(b,ind,orb);
		}

		size_t isThereAnElectronAt(size_t ket1,
		                           size_t ket2,
		                           size_t site,
		                           size_t spin,
		                           size_t orb) const
		{
			if (spin==ProgramGlobals::SPIN_UP)
				return basis1_.isThereAnElectronAt(ket1,site,orb);
			return basis2_.isThereAnElectronAt(ket2,site,orb);
		}

		size_t orbsPerSite(size_t i) const { return orbsPerSite_[i]; }

		size_t orbs() const { return orbsPerSite_[0]; }

		size_t electrons(size_t what) const
		{
			return (what==ProgramGlobals::SPIN_UP) ? basis1_.electrons() : basis2_.electrons();
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

	private:

		bool hasNewPartsCorCdagger(std::pair<size_t,size_t>& newParts,
		                           size_t what,
		                           size_t spin,
		                           const std::pair<size_t,size_t>& orbs) const
		{
			int newPart1=basis1_.electrons();
			int newPart2=basis2_.electrons();

			if (spin==ProgramGlobals::SPIN_UP) newPart1 = basis1_.newPartCorCdagger(what,orbs.first);
			else newPart2 = basis2_.newPartCorCdagger(what,orbs.second);

			if (newPart1<0 || newPart2<0) return false;

			if (newPart1==0 && newPart2==0) return false;
			newParts.first = size_t(newPart1);
			newParts.second = size_t(newPart2);
			return true;
		}

		bool getBra(WordType& bra,
					const WordType& ket1,
					const WordType& ket2,
					size_t what,
					size_t site,
					size_t spin,
					size_t orb) const
		{
			return (spin==ProgramGlobals::SPIN_UP) ? basis1_.getBra(bra,ket1,what,site,orb) :
									 basis2_.getBra(bra,ket2,what,site,orb);
		}

		OrbsPerSite orbsPerSite_;
		BasisType basis1_,basis2_;

	}; // class BasisImmm
} // namespace LanczosPlusPlus
#endif
