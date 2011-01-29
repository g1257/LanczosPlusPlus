
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
#include "BasisOneSpin.h"

namespace LanczosPlusPlus {
	
	class BasisFeAsBasedSc {
	public:
		
		typedef BasisOneSpin BasisType;
		typedef BasisType::WordType WordType;
		enum {SPIN_UP,SPIN_DOWN};
		static size_t const ORBITALS  = BasisType::ORBITALS;
		enum {DESTRUCTOR=BasisType::DESTRUCTOR,CONSTRUCTOR=BasisType::CONSTRUCTOR};
		
		
		BasisFeAsBasedSc(size_t nsite, size_t nup,size_t ndown)
				: basis1_(nsite,nup),basis2_(nsite,ndown)
		{
		}
		
		

		static const WordType& bitmask(size_t i)
		{
			return BasisType::bitmask(i);
		}
		

		size_t size() const { return basis1_.size()*basis2_.size(); }
		

		const WordType& operator()(size_t i,size_t spin) const
		{
			size_t y = i/basis1_.size();
			size_t x = i%basis1_.size();
			return (spin==SPIN_UP) ? basis1_[x] : basis2_[y];
		}
		

		size_t perfectIndex(WordType ket1,WordType ket2) const
		{
			return basis1_.perfectIndex(ket1) + basis2_.perfectIndex(ket2)*basis1_.size();
		}
		

		size_t getN(size_t i,size_t spin,size_t orb) const
		{
			size_t y = i/basis1_.size();
			size_t x = i%basis1_.size();
			return (spin==SPIN_UP) ? basis1_.getN(x,orb) : basis2_.getN(y,orb);
		}
		

		size_t getBraIndex(size_t i,size_t what,size_t sector) const
		{
			size_t y = i/basis1_.size();
			size_t x = i%basis1_.size();
			size_t i1 = basis1_.perfectIndex(basis1_[x]);
			size_t i2 = basis2_.perfectIndex(basis2_[y]);
			size_t spin = sector/2;
			size_t orb = (sector & 1);
			if (spin==SPIN_UP) {
				i1 = basis1_.getBraIndex(x,what,orb);
			} else {
				i2 =  basis2_.getBraIndex(y,what,orb);
			}

			return i1 + i2*basis1_.size();
		}
		

		int doSign(size_t i,size_t site,size_t sector) const
		{
			size_t y = i/basis1_.size();
			size_t x = i%basis1_.size();
			size_t spin = sector/2;
			size_t orb = (sector & 1);
			if (spin==SPIN_UP) {
				return basis1_.doSign(x,site,orb);
			}
			size_t c = basis1_.getN(x);
			int ret = 1;
			if (c&1) ret = -1;
			return ret * basis2_.doSign(y,site,orb);
		}
		
		

	private:
		

		
		
		BasisType basis1_,basis2_;
		
	}; // class BasisFeAsBasedSc
	
} // namespace LanczosPlusPlus
#endif

