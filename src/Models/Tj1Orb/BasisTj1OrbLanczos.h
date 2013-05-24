
/*
*/

#ifndef BASIS_TJ_1ORB_LANCZOS_H
#define BASIS_TJ_1ORB_LANCZOS_H

#include "BitManip.h"
#include "ProgramGlobals.h"

namespace LanczosPlusPlus {
	
	template<typename GeometryType>
	class BasisTj1OrbLanczos {

		typedef ProgramGlobals::PairIntType PairIntType;

	public:		

		typedef ProgramGlobals::WordType WordType;
		typedef PsimagLite::Vector<WordType>::Type VectorWordType;
		static VectorWordType bitmask_;

		enum {OPERATOR_NIL=ProgramGlobals::OPERATOR_NIL,
		      OPERATOR_C=ProgramGlobals::OPERATOR_C,
		      OPERATOR_SZ=ProgramGlobals::OPERATOR_SZ,
		      OPERATOR_CDAGGER=ProgramGlobals::OPERATOR_CDAGGER};

		static int const FERMION_SIGN = -1;

		BasisTj1OrbLanczos(const GeometryType& geometry, SizeType nup,SizeType ndown)
		: geometry_(geometry),nup_(nup),ndown_(ndown)
		{
			assert(bitmask_.size()==0 || bitmask_.size()==geometry_.numberOfSites());
			if (bitmask_.size()==0) doBitmask();
			PsimagLite::Vector<WordType>::Type data1;
			fillOneSector(data1,nup);
			PsimagLite::Vector<WordType>::Type data2;
			fillOneSector(data2,ndown);
			combineAndFilter(data1,data2);
			std::sort(data_.begin(),data_.end());
//			for (SizeType i=0;i<data_.size();i++)
//				std::cout<<"data["<<i<<"]="<<data_[i]<<"\n";
		}
		
		static const WordType& bitmask(SizeType i)
		{
			return bitmask_[i];
		}

		SizeType size() const { return data_.size(); }

		//! Spin up and spin down
		SizeType dofs() const { return 2; }

		SizeType perfectIndex(PsimagLite::Vector<WordType>::Type& kets) const
		{
			assert(kets.size()==2);
			return perfectIndex(kets[0],kets[1]);
		}

		SizeType perfectIndex(WordType ket1,WordType ket2) const
		{
			assert(SizeType(PsimagLite::BitManip::count(ket1))==nup_);
			assert(SizeType(PsimagLite::BitManip::count(ket2))==ndown_);
			SizeType n = geometry_.numberOfSites();
			WordType w = ket2;
			w <<= n;
			w |= ket1;
			SizeType elements = data_.size();
			SizeType i = elements/2;
			SizeType start = 0;
			SizeType end = elements;
			SizeType counter = 0;
			SizeType max = SizeType(0.1*elements);
			if (max>100) max = 100;
		       	if (max<1) max = 1;

			while(counter<max) {
				if (data_[i]==w) return i;
				if (data_[i]>w) {
					if (i<end) end = i;
					i = i/2;
				} else {
					if (i>start) start = i; 
			 		i = (i+elements)/2;
				}
				counter++;
			}
			for (SizeType j=start;j<end;++j) {
				if (data_[j] == w) return j;
			}
			assert(false);
			return 0;
		}

		SizeType electrons(SizeType what) const
		{
			return (what==ProgramGlobals::SPIN_UP) ? nup_ : ndown_;
		}

		WordType operator()(SizeType i,SizeType spin) const
		{
			SizeType n = geometry_.numberOfSites();
			WordType w = data_[i];
			WordType mask = (1<<n);
			mask--;
			if (spin==ProgramGlobals::SPIN_UP) {
				return (w & mask);
			}

			mask <<= n;
			w &= mask;
			w >>= n;
			return w;
		}

		SizeType isThereAnElectronAt(WordType ket1,
		                           WordType ket2,
		                           SizeType site,
		                           SizeType spin) const
		{
			return (spin==ProgramGlobals::SPIN_UP) ? isThereAnElectronAt(ket1,site) : isThereAnElectronAt(ket2,site);
		}

		SizeType getN(WordType ket1,WordType ket2, SizeType site,SizeType spin) const
		{
			return (spin==ProgramGlobals::SPIN_UP) ? getN(ket1,site) : getN(ket2,site);
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
				if (bitmask_[i] & a) s = -s;
				return s;
			}
			int s=(PsimagLite::BitManip::count(a) & 1) ? -1 : 1; // Parity of up
			if (ind==0) return s;

			// ind>0 from now on
			SizeType i = 0;
			SizeType j = ind;
			WordType mask = b;
			mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
			s *= (PsimagLite::BitManip::count(mask) & 1) ? -1 : 1; // Parity of up between i and j
			// Is there a down at i?
			if (bitmask_[i] & b) s = -s;
			return s;
		}

		int doSign(WordType ket1,
		           WordType ket2,
		           SizeType i,
		           SizeType j,
		           SizeType spin) const
		{
			assert(i <= j);
			return (spin==ProgramGlobals::SPIN_UP) ? doSign(ket1,i,j): doSign(ket2,i,j);
		}
		
		int getBra(WordType& bra,
					SizeType operatorLabel,
		            const WordType& ket1,
		            const WordType& ket2,
//		            SizeType what,
		            SizeType site,
		            SizeType spin) const
		{
			switch(operatorLabel) {
			case ProgramGlobals::OPERATOR_C:
			case ProgramGlobals::OPERATOR_CDAGGER:
				return getBraC(bra,ket1,ket2,operatorLabel,site,spin);
			case ProgramGlobals::OPERATOR_SZ:
			case ProgramGlobals::OPERATOR_N:
				return getBraSzOrN(bra,ket1,ket2,operatorLabel,site,spin);
			}
			assert(false);
			return 0;
		}

		PairIntType getBraIndex(const WordType& ket1,
		                        const WordType& ket2,
		                        SizeType operatorLabel,
		                        SizeType site,
		                        SizeType spin,
		                        SizeType orb) const
		{
			WordType bra1 = ket1;
			WordType bra2 = ket2;
			int value = getBra(bra1,operatorLabel,ket1,ket2,site,spin);
			if (value==0) return PairIntType(-1,value);
			if (spin!=ProgramGlobals::SPIN_UP) {
				bra2 = bra1;
				bra1 = ket1;
			}
			int tmp = perfectIndex(bra1,bra2);
			return PairIntType(tmp,value);
		}

		template<typename GeometryType2>
		friend std::ostream& operator<<(std::ostream& os,const BasisTj1OrbLanczos<GeometryType2>& basis);

	private:

		bool isDoublyOccupied(const WordType& ket1,const WordType& ket2) const
		{
			WordType tmp = (ket1 & ket2);
			return (tmp>0);
		}

		void fillOneSector(PsimagLite::Vector<WordType>::Type& data1,SizeType npart) const
		{
			/* compute size of basis */
			SizeType hilbert=1;
			int n=geometry_.numberOfSites();
			SizeType m=1;
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
			for (SizeType i=0;i<hilbert;i++) {
				data1[i] = ket;
				n=m=0;
				for (;(ket&3)!=1;n++,ket>>=1) {
					m += ket&1;
				}
				ket = ((ket+1)<<n) ^ ((1<<m)-1);
			}
		}

		void combineAndFilter(const PsimagLite::Vector<WordType>::Type& data1,const PsimagLite::Vector<WordType>::Type& data2)
		{
			WordType tmp = 0;
			WordType tmp2 = 0;
			SizeType n=geometry_.numberOfSites();
			for (SizeType i=0;i<data1.size();i++) {
				for (SizeType j=0;j<data2.size();j++) {
					tmp = (data1[i] & data2[j]);
					if (tmp>0) continue; // there's one or more double occupied
					tmp2 = data2[j];
					tmp2 <<= n;
					data_.push_back(tmp2 | data1[i]);
				}
			}
		}

		SizeType isThereAnElectronAt(WordType ket,SizeType site) const
		{
			return (ket & bitmask_[site]) ? 1 : 0;
		}

		SizeType getN(WordType ket,SizeType site) const
		{
			return isThereAnElectronAt(ket,site);
		}

		void doBitmask()
		{
			SizeType n = geometry_.numberOfSites();
			bitmask_.resize(n);
			bitmask_[0]=1ul;
			for (SizeType i=1;i<n;i++)
				bitmask_[i] = bitmask_[i-1]<<1;
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

		int getBraC(WordType& bra,
					const WordType& ket1,
					const WordType& ket2,
					SizeType what,
					SizeType site,
					SizeType spin) const
		{
			if (spin==ProgramGlobals::SPIN_UP) {
				int b1 = getBraC(bra,ket1,what,site);
				if (b1==0) return 0;
				return (isDoublyOccupied(bra,ket2)) ? 0 : 1;
			}

			int b2 = getBraC(bra,ket2,what,site);
			if (b2==0) return 0;
			return (isDoublyOccupied(ket1,bra)) ? 0 : 1;
		}

		int getBraC(WordType& bra,const WordType& ket,SizeType what,SizeType site) const
		{
			WordType si=(ket & bitmask_[site]);
			if (what==OPERATOR_C) {
				if (si>0) {
					bra = (ket ^ bitmask_[site]);
				} else {
					return 0; // cannot destroy, there's nothing
				}
			} else {
				if (si==0) {
					bra = (ket ^ bitmask_[site]);
				} else {
					return 0; // cannot construct, there's already one
				}
			}
			return 1;
		}

		int getBraSzOrN(WordType& bra,
		                const WordType& ket1,
		                const WordType& ket2,
		                SizeType operatorLabel,
		                SizeType site,
		                SizeType spin) const
		{
			assert(spin==ProgramGlobals::SPIN_UP); // spin index is bogus here
			if (spin==ProgramGlobals::SPIN_UP) bra = ket1;
			else bra = ket2;
			WordType siup=(ket1 & bitmask_[site]);
			WordType sidown=(ket2 & bitmask_[site]);
			if (siup>0) siup=1;
			if (sidown>0) sidown=1;
			int isign = (operatorLabel == OPERATOR_SZ) ? -1 : 1;
			return (siup + isign * sidown);
		}

		const GeometryType& geometry_;
		SizeType nup_;
		SizeType ndown_;
		PsimagLite::Vector<WordType>::Type data_;
	}; // class BasisTj1OrbLanczos
	
	template<typename GeometryType>
	std::ostream& operator<<(std::ostream& os,const BasisTj1OrbLanczos<GeometryType>& basis)
	{
		for (SizeType i=0;i<basis.data_.size();i++)
			os<<i<<" "<<basis.data_[i]<<"\n";
		return os;
	}

	template<typename GeometryType>
	typename BasisTj1OrbLanczos<GeometryType>::VectorWordType
	BasisTj1OrbLanczos<GeometryType>::bitmask_;


} // namespace LanczosPlusPlus
#endif

