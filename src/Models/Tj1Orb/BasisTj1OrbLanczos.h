
/*
*/

#ifndef BASIS_TJ_1ORB_LANCZOS_H
#define BASIS_TJ_1ORB_LANCZOS_H

#include "BitManip.h"

namespace LanczosPlusPlus {
	
	template<typename GeometryType>
	class BasisTj1OrbLanczos {

	public:		

		typedef unsigned int long long WordType;

		static std::vector<WordType> bitmask_;

		enum {SPIN_UP,SPIN_DOWN};

		enum {DESTRUCTOR,CONSTRUCTOR};

		static int const FERMION_SIGN = -1;

		BasisTj1OrbLanczos(const GeometryType& geometry, size_t nup,size_t ndown)
		: geometry_(geometry),nup_(nup),ndown_(ndown)
		{
			doBitmask();
			std::vector<WordType> data1;
			fillOneSector(data1,nup);
			std::vector<WordType> data2;
			fillOneSector(data2,ndown);
			combineAndFilter(data1,data2);
//			for (size_t i=0;i<data_.size();i++)
//				std::cout<<"data["<<i<<"]="<<data_[i]<<"\n";
		}
		
		static const WordType& bitmask(size_t i)
		{
			return bitmask_[i];
		}

		size_t size() const { return data_.size(); }

		//! Spin up and spin down
		size_t dofs() const { return 2; }

		size_t perfectIndex(std::vector<WordType>& kets) const
		{
			assert(kets.size()==2);
			return perfectIndex(kets[0],kets[1]);
		}

		size_t perfectIndex(WordType ket1,WordType ket2) const
		{
			size_t n = geometry_.numberOfSites();
			WordType w = ket2;
			w <<= n;
			w |= ket1;
			for (size_t i=0;i<data_.size();i++) {
				if (data_[i] == w) return i;
			}
			assert(false);
			return 0;
		}

		size_t electrons(size_t what) const
		{
			return (what==SPIN_UP) ? nup_ : ndown_;
		}

		WordType operator()(size_t i,size_t spin) const
		{
			size_t n = geometry_.numberOfSites();
			WordType w = data_[i];
			WordType mask = (1<<n);
			mask--;
			if (spin==SPIN_UP) {
				return (w & mask);
			}

			mask <<= n;
			w &= mask;
			w >>= n;
			return w;
		}

		size_t isThereAnElectronAt(WordType ket1,
		                           WordType ket2,
		                           size_t site,
		                           size_t spin) const
		{
			return (spin==SPIN_UP) ? isThereAnElectronAt(ket1,site) : isThereAnElectronAt(ket2,site);
		}

		size_t getN(WordType ket1,WordType ket2, size_t site,size_t spin) const
		{
			return (spin==SPIN_UP) ? getN(ket1,site) : getN(ket2,site);
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
				if (bitmask_[i] & a) s = -s;
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
			if (bitmask_[i] & b) s = -s;
			return s;
		}

		int doSign(WordType ket1,
		           WordType ket2,
		           size_t i,
		           size_t j,
		           size_t spin) const
		{
			assert(i <= j);
			return (spin==SPIN_UP) ? doSign(ket1,i,j): doSign(ket2,i,j);
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
			return (spin==SPIN_UP) ? getBra(bra,ket1,what,site) :
									 getBra(bra,ket2,what,site);
		}

	private:

		void fillOneSector(std::vector<WordType>& data1,size_t npart) const
		{
			/* compute size of basis */
			size_t hilbert=1;
			int n=geometry_.numberOfSites();
			size_t m=1;
			for (;m<=npart;n--,m++)
				hilbert=hilbert*n/m;

			if (data1.size()!=hilbert) {
				data1.clear();
				data1.resize(hilbert);
			}

			if (npart==0) {
				data1[0]=0;
				return;
			}

			/* define basis states */
			WordType ket = (1ul<<npart)-1;
			for (size_t i=0;i<hilbert;i++) {
				data1[i] = ket;
				n=m=0;
				for (;(ket&3)!=1;n++,ket>>=1) {
					m += ket&1;
				}
				ket = ((ket+1)<<n) ^ ((1<<m)-1);
			}
		}

		void combineAndFilter(const std::vector<WordType>& data1,const std::vector<WordType>& data2)
		{
			WordType tmp = 0;
			WordType tmp2 = 0;
			size_t n=geometry_.numberOfSites();
			for (size_t i=0;i<data1.size();i++) {
				for (size_t j=0;j<data2.size();j++) {
					tmp = (data1[i] & data2[j]);
					if (tmp>0) continue; // there's one or more double occupied
					tmp2 = data2[j];
					tmp2 <<= n;
					data_.push_back(tmp2 | data1[i]);
				}
			}
		}

		size_t isThereAnElectronAt(WordType ket,size_t site) const
		{
			return (ket & bitmask_[site]) ? 1 : 0;
		}

		size_t getN(WordType ket,size_t site) const
		{
			return isThereAnElectronAt(ket,site);
		}

		void doBitmask()
		{
			size_t n = geometry_.numberOfSites();
			bitmask_.resize(n);
			bitmask_[0]=1ul;
			for (size_t i=1;i<n;i++)
				bitmask_[i] = bitmask_[i-1]<<1;
		}

		int doSign(WordType ket,size_t i,size_t j) const
		{
			assert(i <= j);
			size_t x0 = (i+1); // i+1 cannot be the last site, 'cause i<j
			size_t x1 = j;

			size_t sum = getNbyKet(ket,x0,x1);

			// at site i we need to be carefull
			x0 = i;
			x1 = (i+1);
			sum += getNbyKet(ket,x0,x1);

			// same at site j
			x0 = j;
			x1 = j;
			sum += getNbyKet(ket,x0,x1);

			return (sum & 1) ? FERMION_SIGN : 1;
		}

		size_t getNbyKet(size_t ket,size_t from,size_t upto) const
		{
			size_t sum = 0;
			size_t counter = from;
			while(counter<upto) {
				if (ket & bitmask_[counter]) sum++;
				counter++;
			}
			return sum;
		}

		bool getBra(WordType& bra, const WordType& ket,size_t what,size_t site) const
		{

			WordType si=(ket & bitmask_[site]);
			if (what==DESTRUCTOR) {
				if (si>0) {
					bra = (ket ^ bitmask_[site]);
				} else {
					return false; // cannot destroy, there's nothing
				}
			} else {
				if (si==0) {
					bra = (ket ^ bitmask_[site]);
				} else {
					return false; // cannot construct, there's already one
				}
			}
			return true;
		}

		const GeometryType& geometry_;
		size_t nup_;
		size_t ndown_;
		std::vector<WordType> data_;
	}; // class BasisTj1OrbLanczos
	
	template<typename GeometryType>
	std::vector<typename BasisTj1OrbLanczos<GeometryType>::WordType> BasisTj1OrbLanczos<GeometryType>::bitmask_;

} // namespace LanczosPlusPlus
#endif

