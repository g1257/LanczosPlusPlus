
/*
*/

#ifndef BASISHUBBARDLANCZOS_H
#define BASISHUBBARDLANCZOS_H
#include "BasisOneSpin.h"

namespace LanczosPlusPlus {
	
	template<typename GeometryType>
	class BasisHubbardLanczos {
	public:
		
		typedef BasisOneSpin BasisType;
		typedef typename BasisType::WordType WordType;

		enum {SPIN_UP,SPIN_DOWN};

		enum {DESTRUCTOR=BasisType::DESTRUCTOR,CONSTRUCTOR=BasisType::CONSTRUCTOR};

		static int const FERMION_SIGN = BasisType::FERMION_SIGN;

		BasisHubbardLanczos(const GeometryType& geometry, size_t nup,size_t ndown) 
		: basis1_(geometry.numberOfSites(),nup),
		  basis2_(geometry.numberOfSites(),ndown)
		{} 
		
		static const WordType& bitmask(size_t i)
		{
			return BasisType::bitmask(i);
		}

		size_t size() const { return basis1_.size()*basis2_.size(); }

		size_t perfectIndex(WordType ket1,WordType ket2) const
		{
			return basis1_.perfectIndex(ket1) + 
			       basis2_.perfectIndex(ket2)*basis1_.size();
		}

		size_t electrons(size_t what) const
		{
			return (what==SPIN_UP) ? basis1_.electrons() : basis2_.electrons();
		}

		const WordType& operator()(size_t i,size_t spin) const
		{
			size_t y = i/basis1_.size();
			size_t x = i%basis1_.size();
			return (spin==SPIN_UP) ? basis1_[x] : basis2_[y];
		}

		size_t isThereAnElectronAt(WordType ket1,
		                           WordType ket2,
		                           size_t site,
		                           size_t spin) const
		{
			if (spin==SPIN_UP)
				return basis1_.isThereAnElectronAt(ket1,site);
			return basis2_.isThereAnElectronAt(ket2,site);
		}

		size_t getN(WordType ket1,WordType ket2, size_t site,size_t spin) const
		{
			return (spin==SPIN_UP) ? basis1_.getN(ket1,site) : basis2_.getN(ket2,site);
		}
		
		int doSignGf(WordType a, WordType b,size_t ind,size_t sector) const
		{
			if (sector==SPIN_UP) {
				if (ind==0) return 1;

				// ind>0 from now on
				size_t i = 0;
				size_t j = ind;
				WordType mask = a;
				mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
				int s=(PsimagLite::BitManip::count(mask) & 1) ? -1 : 1; // Parity of up between i and j
				// Is there an up at i?
				if (BasisType::bitmask(i) & a) s = -s;
				return s;
			}
			int s=(PsimagLite::BitManip::count(a) & 1) ? -1 : 1; // Parity of up
			if (ind==0) return s;

			// ind>0 from now on
			size_t i = 0;
			size_t j = ind;
			WordType mask = b;
			mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
			s=(PsimagLite::BitManip::count(mask) & 1) ? -1 : 1; // Parity of up between i and j
			// Is there a down at i?
			if (BasisType::bitmask(i) & b) s = -s;
			return s;
		}

		int doSign(WordType ket1,
		           WordType ket2,
		           size_t i,
		           size_t j,
		           size_t spin) const
		{
			assert(i <= j);
			return (spin==SPIN_UP) ? basis1_.doSign(ket1,i,j): basis2_.doSign(ket2,i,j);
		}
		
		int doSign(WordType ket,
		           size_t i,
		           size_t spin) const
		{
			return (spin==SPIN_UP) ? basis1_.doSign(ket,i): basis2_.doSign(ket,i);
		}
		
		int getBraIndex(const WordType& ket1,
		                   const WordType& ket2,
		                   size_t what,
		                   size_t site,
		                   size_t spin) const
		{
			WordType bra;
			bool b = getBra(bra,ket1,ket2,what,site,spin);
			if (!b) return -1;
			return (spin==SPIN_UP) ? perfectIndex(bra,ket2) :
			                         perfectIndex(ket1,bra);
		}
		
		bool getBra(WordType& bra,
		            const WordType& ket1,
		            const WordType& ket2,
		            size_t what,
		            size_t site,
		            size_t spin) const
		{
			return (spin==SPIN_UP) ? basis1_.getBra(bra,ket1,what,site) :
			                         basis2_.getBra(bra,ket2,what,site);
		}

	private:

		BasisType basis1_,basis2_;
		
	}; // class BasisHubbardLanczos
	
} // namespace LanczosPlusPlus
#endif

