
/*
*/

#ifndef BASISHUBBARDLANCZOS_H
#define BASISHUBBARDLANCZOS_H
#include "Matrix.h"

namespace LanczosPlusPlus {
	
	class BasisHubbardLanczos {
	public:
		
		typedef unsigned int long long WordType;
		static size_t nsite_;
		static PsimagLite::Matrix<size_t> comb_;
		static std::vector<WordType> bitmask_; 
		
		BasisHubbardLanczos(size_t nsite, size_t npart) : npart_(npart)
		{
			if (nsite_>0 && nsite!=nsite_)
				throw std::runtime_error("BasisHubbardLanczos: All basis must have same number of sites\n");
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
		

	private:
		

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
		
	}; // class BasisHubbardLanczos

	size_t BasisHubbardLanczos::nsite_=0;
	PsimagLite::Matrix<size_t> BasisHubbardLanczos::comb_;
	std::vector<typename BasisHubbardLanczos::WordType> BasisHubbardLanczos::bitmask_; 
	
} // namespace LanczosPlusPlus
#endif

