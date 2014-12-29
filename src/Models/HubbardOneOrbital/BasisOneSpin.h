
/*
*/

#ifndef BASIS_ONE_SPIN_H
#define BASIS_ONE_SPIN_H
#include "Matrix.h"
#include "BitManip.h"
#include "ProgramGlobals.h"

namespace LanczosPlusPlus {

	class BasisOneSpin {

	public:

		static int const FERMION_SIGN  = -1;
		typedef ProgramGlobals::WordType WordType;
		static SizeType nsite_;
		static PsimagLite::Matrix<SizeType> comb_;
		static PsimagLite::Vector<WordType>::Type bitmask_;

		enum {OPERATOR_NIL=ProgramGlobals::OPERATOR_NIL,
		      OPERATOR_C=ProgramGlobals::OPERATOR_C,
		      OPERATOR_SZ=ProgramGlobals::OPERATOR_SZ,
		      OPERATOR_CDAGGER=ProgramGlobals::OPERATOR_CDAGGER};

		BasisOneSpin(SizeType nsite, SizeType npart)
		: npart_(npart)
		{
			if (nsite_>0 && nsite!=nsite_)
				throw std::runtime_error("BasisOneSpin: All basis must have same number of sites\n");
			nsite_ = nsite;
			doCombinatorial();
			doBitmask();

			/* compute size of basis */
			SizeType hilbert=1;
			int n=nsite;
			SizeType m=1;
			for (;m<=npart;n--,m++)
				hilbert=hilbert*n/m;

			if (data_.size()!=hilbert) {
				data_.clear();
				data_.resize(hilbert);
			}

			if (npart==0) {
				data_[0]=0;
				size_ = 1;
				return;
			}

			/* define basis states */
			WordType ket = (1ul<<npart)-1;
			for (SizeType i=0;i<hilbert;i++) {
				data_[i] = ket;
				n=m=0;
				for (;(ket&3)!=1;n++,ket>>=1) {
					m += ket&1;
				}
				ket = ((ket+1)<<n) ^ ((1<<m)-1);
			}
			size_ = hilbert;
		}


		SizeType size() const { return size_; }

		const WordType& operator[](SizeType i) const
		{
			return data_[i];
		}

		SizeType perfectIndex(WordType state) const
		{
			SizeType n=0;
			for (SizeType b=0,c=1;state>0;b++,state>>=1)
				if (state&1) n += comb_(b,c++);

			assert(n<data_.size());
			return n;
		}

		static const WordType& bitmask(SizeType i)
		{
			return bitmask_[i];
		}

		SizeType electrons() const { return npart_; }

		SizeType isThereAnElectronAt(WordType ket,SizeType site) const
		{
			return (ket & bitmask_[site]) ? 1 : 0;
		}

		SizeType getN(WordType ket,SizeType site) const
		{
			return isThereAnElectronAt(ket,site);
		}

		int doSign(WordType a, SizeType i) const
		{
			if (i==nsite_-1) return 1;

			a &= ((1 << (i+1)) - 1) ^ ((1 << nsite_) - 1);
			// Parity of single occupied between i and nsite-1
			int s=(PsimagLite::BitManip::count(a) & 1) ? FERMION_SIGN : 1;
			return s;
		}

		int doSign(WordType ket,SizeType i,SizeType j) const
		{
			assert(i <= j);
			SizeType x0 = (i+1); // i+1 cannot be the last site, 'cause i<j
			SizeType x1 = j;

			SizeType sum = getNbyKet(ket,x0,x1);

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

		bool getBra(WordType& bra, const WordType& ket,SizeType what,SizeType site) const
		{
			WordType si=(ket & bitmask_[site]);
			if (what==OPERATOR_C) {
				if (si>0) {
					bra = (ket ^ bitmask_[site]);
					return true;
				} else {
					return false; // cannot destroy, there's nothing
				}
			} else if (what==OPERATOR_CDAGGER) {
				if (si==0) {
					bra = (ket ^ bitmask_[site]);
					return true;
				} else {
					return false; // cannot construct, there's already one
				}
			} else if (what==ProgramGlobals::OPERATOR_N) {
				if (si==0) return false;
				bra = ket;
				return true;
			}
			PsimagLite::String str = ProgramGlobals::unknownOperator(what);
			throw std::runtime_error(str.c_str());
		}

		void print(std::ostream& os) const
		{
			SizeType hilbert = 1;
			hilbert <<= nsite_;
			ProgramGlobals::printBasisVector(os,hilbert,data_);
		}

	private:

		SizeType getNbyKet(SizeType ket,SizeType from,SizeType upto) const
		{
			SizeType sum = 0;
			SizeType counter = from;
			while(counter<upto) {
				if (ket & bitmask_[counter]) sum++;
				counter++;
			}
			return sum;
		}

// 		SizeType getNbyKet(SizeType ket) const
// 		{
// 			SizeType sum = 0;
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

			for (SizeType n=0;n<nsite_;n++)
				for (SizeType i=0;i<nsite_;i++)
					comb_(n,i)=0;

			for (SizeType n=0;n<nsite_;n++) {
				SizeType m = 0;
				int j = n;
				SizeType i = 1;
				SizeType cnm  = 1;
				for (;m<=n/2;m++,cnm=cnm*j/i,i++,j--)
					comb_(n,m) = comb_(n,n-m) = cnm;
			}
		}

		void doBitmask()
		{
			bitmask_.resize(nsite_);
			bitmask_[0]=1ul;
			for (SizeType i=1;i<nsite_;i++)
				bitmask_[i] = bitmask_[i-1]<<1;
		}


		SizeType size_;
		SizeType npart_;
		PsimagLite::Vector<WordType>::Type data_;

	}; // class BasisOneSpin

	SizeType BasisOneSpin::nsite_=0;

	PsimagLite::Matrix<SizeType> BasisOneSpin::comb_;

	PsimagLite::Vector<BasisOneSpin::WordType>::Type BasisOneSpin::bitmask_;

} // namespace LanczosPlusPlus
#endif // BASIS_ONE_SPIN_H

