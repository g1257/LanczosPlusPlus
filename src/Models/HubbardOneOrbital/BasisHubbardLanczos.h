
/*
*/

#ifndef BASISHUBBARDLANCZOS_H
#define BASISHUBBARDLANCZOS_H
#include "BasisOneSpin.h"
#include "BasisBase.h"

namespace LanczosPlusPlus {
	
	template<typename GeometryType>
	class BasisHubbardLanczos : public BasisBase<GeometryType> {

		typedef ProgramGlobals::PairIntType PairIntType;

	public:
		
		typedef BasisOneSpin BasisType;
		typedef BasisBase<GeometryType> BaseType;
		typedef typename BaseType::WordType WordType;
		typedef typename BaseType::VectorWordType VectorWordType;

		static int const FERMION_SIGN = BasisType::FERMION_SIGN;

		BasisHubbardLanczos(const GeometryType& geometry, SizeType nup,SizeType ndown) 
		: basis1_(geometry.numberOfSites(),nup),
		  basis2_(geometry.numberOfSites(),ndown)
		{} 
		
		static const WordType& bitmask(SizeType i)
		{
			return BasisType::bitmask(i);
		}

		SizeType size() const { return basis1_.size()*basis2_.size(); }

		//! Spin up and spin down
		SizeType dofs() const { return 2; }

		SizeType perfectIndex(const VectorWordType& kets) const
		{
			assert(kets.size()==2);
			return perfectIndex(kets[0],kets[1]);
		}

		SizeType perfectIndex(WordType ket1,WordType ket2) const
		{
			return basis1_.perfectIndex(ket1) + 
			       basis2_.perfectIndex(ket2)*basis1_.size();
		}

		virtual SizeType perfectIndex(WordType newKet,
	                                      SizeType ispace,
	                                      SizeType spinOfNew) const
		{
			throw PsimagLite::RuntimeError("perfectIndex\n");
		}

		SizeType electrons(SizeType what) const
		{
			return (what==ProgramGlobals::SPIN_UP) ? basis1_.electrons() : basis2_.electrons();
		}

		WordType operator()(SizeType i,SizeType spin) const
		{
			SizeType y = i/basis1_.size();
			SizeType x = i%basis1_.size();
			return (spin==ProgramGlobals::SPIN_UP) ? basis1_[x] : basis2_[y];
		}

		SizeType isThereAnElectronAt(WordType ket1,
		                             WordType ket2,
		                             SizeType site,
		                             SizeType spin,
		                             SizeType orb) const
		{
			if (spin==ProgramGlobals::SPIN_UP)
				return basis1_.isThereAnElectronAt(ket1,site);
			return basis2_.isThereAnElectronAt(ket2,site);
		}

		SizeType getN(WordType ket1,WordType ket2, SizeType site,SizeType spin, SizeType orb) const
		{
			return (spin==ProgramGlobals::SPIN_UP) ? basis1_.getN(ket1,site) : basis2_.getN(ket2,site);
		}
		
		int doSignGf(WordType a, WordType b,SizeType ind,SizeType sector,SizeType orb) const
		{
			if (sector==ProgramGlobals::SPIN_UP) {
				if (ind==0) return 1;

				// ind>0 from now on
				SizeType i = 0;
				SizeType j = ind;
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
			SizeType i = 0;
			SizeType j = ind;
			WordType mask = b;
			mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
			s=(PsimagLite::BitManip::count(mask) & 1) ? -1 : 1; // Parity of up between i and j
			// Is there a down at i?
			if (BasisType::bitmask(i) & b) s = -s;
			return s;
		}

		int doSign(WordType ket1,
		           WordType ket2,
		           SizeType i,
		           SizeType orb1,
		           SizeType j,
		           SizeType orb2,
		           SizeType spin) const
		{
			assert(i <= j);
			return (spin==ProgramGlobals::SPIN_UP) ? basis1_.doSign(ket1,i,j): basis2_.doSign(ket2,i,j);
		}

		PairIntType getBraIndex(WordType ket1,
		                        WordType ket2,
		                        SizeType what,
		                        SizeType site,
		                        SizeType spin,
		                        SizeType orb) const
		{
			WordType bra = 0;
			bool b = getBra(bra,ket1,ket2,what,site,spin);
			if (!b) return PairIntType(-1,1);
			int tmp = (spin==ProgramGlobals::SPIN_UP) ? perfectIndex(bra,ket2) :
			                         perfectIndex(ket1,bra);
			return PairIntType(tmp,1);
		}
		
		bool getBra(WordType& bra,
		            const WordType& ket1,
		            const WordType& ket2,
		            SizeType what,
		            SizeType site,
		            SizeType spin) const
		{
			return (spin==ProgramGlobals::SPIN_UP) ? basis1_.getBra(bra,ket1,what,site) :
									 basis2_.getBra(bra,ket2,what,site);
		}

		SizeType orbsPerSite(SizeType i) const { return 1; } 
		
		SizeType orbs() const { return 1; }

	private:

		BasisType basis1_,basis2_;
		
	}; // class BasisHubbardLanczos
	
} // namespace LanczosPlusPlus
#endif

