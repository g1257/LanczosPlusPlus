
/*
*/

#ifndef BASIS_ONE_SPIN_H
#define BASIS_ONE_SPIN_H
#include "Matrix.h"
#include "BitManip.h"

namespace LanczosPlusPlus {
	
	class BasisOneSpin {

	public:
		
		static int const FERMION_SIGN  = -1;
		typedef unsigned int long long WordType;
		static size_t nsite_;
		static PsimagLite::Matrix<size_t> comb_;
		static std::vector<WordType> bitmask_; 

		enum {DESTRUCTOR,CONSTRUCTOR};

		BasisOneSpin(size_t nsite, size_t npart) 
		: npart_(npart)
		{
			if (nsite_>0 && nsite!=nsite_)
				throw std::runtime_error("BasisOneSpin: All basis must have same number of sites\n");
			nsite_ = nsite;
			doCombinatorial();
			doBitmask();

			/* compute size of basis */
			size_t hilbert=1;
			int n=nsite;
			size_t m=1;
			for (;m<=npart;n--,m++)
				hilbert=hilbert*n/m;

			if (data_.size()!=hilbert) {
				data_.clear();
				data_.resize(hilbert);
			}

			if (npart==0) {
				data_[0]=0;
				return;
			}
			
			/* define basis states */
			WordType ket = (1ul<<npart)-1;
			for (size_t i=0;i<hilbert;i++) {
				data_[i] = ket;
				n=m=0;
				for (;(ket&3)!=1;n++,ket>>=1) {
					m += ket&1;
				}
				ket = ((ket+1)<<n) ^ ((1<<m)-1);
			}
			size_ = hilbert;
		} 
		

		size_t size() const { return size_; } 

		const WordType& operator[](size_t i) const
		{
			return data_[i];
		} 

		size_t perfectIndex(WordType state) const
		{
			size_t n=0;
			for (size_t b=0,c=1;state>0;b++,state>>=1)
				if (state&1) n += comb_(b,c++);

			return n;
		} 

		static const WordType& bitmask(size_t i)
		{
			return bitmask_[i];
		}

		size_t electrons() const { return npart_; }

		size_t isThereAnElectronAt(size_t ket,size_t site) const
		{
			return (ket & bitmask_[site]) ? 1 : 0;
		}
		
		size_t getN(size_t i) const
		{
			return isThereAnElectronAt(data_[i],i);
		}

		int doSign(WordType a, size_t i) const
		{
			if (i==nsite_-1) return 1;

			a &= ((1 << (i+1)) - 1) ^ ((1 << nsite_) - 1);
			// Parity of single occupied between i and nsite-1
			int s=(PsimagLite::BitManip::count(a) & 1) ? FERMION_SIGN : 1;
			return s;
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
					return false; // cannot contruct, there's already one
				}
			}
			return true;
		}
	
	private:

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

// 		size_t getNbyKet(size_t ket) const
// 		{
// 			size_t sum = 0;
// 			WordType ketCopy = ket;
// 			while(ketCopy) {
// 				if (ketCopy & 1) sum++;
// 				ketCopy <<= 1;
// 			}
// 			return sum;
// 		}

		void doCombinatorial()
		{
			/* look-up table for binomial coefficients */
			comb_.reset(nsite_,nsite_);

			for (size_t n=0;n<nsite_;n++)
				for (size_t i=0;i<nsite_;i++)
					comb_(n,i)=0;

			for (size_t n=0;n<nsite_;n++) {
				size_t m = 0;
				int j = n;
				size_t i = 1;
				size_t cnm  = 1;
				for (;m<=n/2;m++,cnm=cnm*j/i,i++,j--)
					comb_(n,m) = comb_(n,n-m) = cnm;
			}
		} 

		void doBitmask()
		{
			bitmask_.resize(nsite_);
			bitmask_[0]=1ul;
			for (size_t i=1;i<nsite_;i++)
				bitmask_[i] = bitmask_[i-1]<<1;
		} 
		
		
		size_t size_;
		size_t npart_;
		std::vector<WordType> data_;
		
	}; // class BasisOneSpin

	size_t BasisOneSpin::nsite_=0;

	PsimagLite::Matrix<size_t> BasisOneSpin::comb_;

	std::vector<typename BasisOneSpin::WordType> BasisOneSpin::bitmask_; 

} // namespace LanczosPlusPlus
#endif // BASIS_ONE_SPIN_H

