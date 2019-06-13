
/*
*/

#ifndef BASIS_ONE_SPIN_H
#define BASIS_ONE_SPIN_H
#include "Matrix.h"
#include "BitManip.h"
#include "ProgramGlobals.h"
#include "LabeledOperator.h"

namespace LanczosPlusPlus {

	class BasisOneSpin {

	public:

		typedef ProgramGlobals::WordType WordType;
		typedef LabeledOperator LabeledOperatorType;

		static int const FERMION_SIGN  = -1;
		static SizeType nsite_;
		static PsimagLite::Matrix<SizeType> comb_;

		BasisOneSpin(SizeType nsite, SizeType npart)
		: npart_(npart)
		{
			if (nsite_>0 && nsite!=nsite_)
				throw std::runtime_error("BasisOneSpin: All basis must have same number of sites\n");
			nsite_ = nsite;
			doCombinatorial();
			ProgramGlobals::doBitmask(nsite_);

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
			return ProgramGlobals::bitmask(i);
		}

		SizeType electrons() const { return npart_; }

		static SizeType isThereAnElectronAt(WordType ket,SizeType site)
		{
			return (ket & ProgramGlobals::bitmask(site)) ? 1 : 0;
		}

		static SizeType getN(WordType ket,SizeType site)
		{
			return isThereAnElectronAt(ket,site);
		}

		static int doSign(WordType ket,SizeType i,SizeType j)
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

		bool getBra(WordType& bra,
		            const WordType& ket,
		            const LabeledOperatorType& lOperator,
		            SizeType site) const
		{
			WordType si=(ket & ProgramGlobals::bitmask(site));
			if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_C) {
				if (si>0) {
					bra = (ket ^ ProgramGlobals::bitmask(site));
					return true;
				} else {
					return false; // cannot destroy, there's nothing
				}
			} else if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_CDAGGER) {
				if (si==0) {
					bra = (ket ^ ProgramGlobals::bitmask(site));
					return true;
				} else {
					return false; // cannot construct, there's already one
				}
			} else if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_N) {
				if (si==0) return false;
				bra = ket;
				return true;
			}

			PsimagLite::String str = lOperator.unknownOperator();
			throw std::runtime_error(str.c_str());
		}

		void print(std::ostream& os,bool isBinary) const
		{
			SizeType hilbert = 1;
			hilbert <<= nsite_;
			if (isBinary) {
				ProgramGlobals::printBasisBinary(os,hilbert,data_);
			} else {
				ProgramGlobals::printBasisDecimal(os,40,data_);
			}
		}

		static SizeType comb(SizeType n, SizeType m) { return comb_(n, m); }

	private:

		static SizeType getNbyKet(SizeType ket,SizeType from,SizeType upto)
		{
			SizeType sum = 0;
			SizeType counter = from;
			while(counter<upto) {
				if (ket & ProgramGlobals::bitmask(counter)) sum++;
				counter++;
			}

			return sum;
		}

		static void doCombinatorial()
		{
			/* look-up table for binomial coefficients */
			comb_.resize(2*nsite_ + 2, 2*nsite_ + 2, 0);

			for (SizeType n=0;n<comb_.n_row();n++) {
				SizeType m = 0;
				int j = n;
				SizeType i = 1;
				SizeType cnm  = 1;
				for (;m<=n/2;m++,cnm=cnm*j/i,i++,j--)
					comb_(n,m) = comb_(n,n-m) = cnm;
			}
		}

		SizeType size_;
		SizeType npart_;
		PsimagLite::Vector<WordType>::Type data_;

	}; // class BasisOneSpin
} // namespace LanczosPlusPlus
#endif // BASIS_ONE_SPIN_H

