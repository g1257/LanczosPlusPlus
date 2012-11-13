
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

#ifndef BASIS_ONE_SPIN_FE_AS_H
#define BASIS_ONE_SPIN_FE_AS_H
#include "Matrix.h"
#include "BitManip.h"
#include "Partitions.h"

namespace LanczosPlusPlus {
	
	class BasisOneSpinFeAs {

		typedef Partitions PartitionsType;

	public:
		
		typedef unsigned int long long WordType;
		enum {DESTRUCTOR,CONSTRUCTOR};
		static size_t orbitals_;
		static int const FERMION_SIGN  = -1;
		static size_t nsite_;
		static PsimagLite::Matrix<size_t> comb_;
		static std::vector<WordType> bitmask_; 
		
		BasisOneSpinFeAs(size_t nsite, size_t npart,size_t orbitals)
				: npart_(npart)
		{
			if (nsite_>0 && nsite!=nsite_)
				throw std::runtime_error("BasisOneSpinFeAs: All basis must have same number of sites\n");
			orbitals_=orbitals;
			nsite_ = nsite;
			doCombinatorial();
			doBitmask();

			/* compute size of basis */
			if (npart==0) {
				data_.resize(1);
				size_ =1;
				data_[0]=0;
				return;
			}
			size_ = 0;
			for (size_t na=0;na<=npart;na++) {
					size_t nb = npart - na;
					size_ += comb_(nsite_,na) * comb_(nsite_,nb);
			}
			data_.resize(size_);

			// compute basis:
			size_t counter = 0;
			PartitionsType partitions(npart,orbitals_);
			for (size_t i=0;i<partitions.size();i++) {
				const std::vector<size_t>& na = partitions(i);
				std::vector<std::vector<WordType> > basisA(orbitals_);
				for (size_t orb=0;orb<orbitals_;orb++) {
					fillPartialBasis(basisA[orb],na[orb]);
				}
				collateBasis(counter,basisA);
			}
		}

		size_t size() const { return size_; }

		const WordType& operator[](size_t i) const
		{
			return data_[i];
		}

		size_t perfectIndex(WordType ket) const
		{
			for (size_t i=0;i<data_.size();i++)
				if (data_[i]==ket) return i;
			throw std::runtime_error("perfectindex\n");
		}

//			WordType ketA=0,ketB=0;
//			uncollateKet(ketA,ketB,ket);
//			// p(ket) = \sum_{na'=0}^{na'<na} S_na' * S_nb'
//			//			+ p_A(ket_A)*S_nb + p_B(ket_B)
//			// where S_x = C^n_x
//			size_t na = PsimagLite::BitManip::count(ketA);
//			// note nb = PsimagLite::BitManip::count(ketB)
//			// or nb  = npart -na
//			size_t s = 0;
//			for (size_t nap=0;nap<na;nap++) {
//				size_t nbp = npart_ - nap;
//				s += comb_(nsite_,nap) * comb_(nsite_,nbp);
//			}
//			size_t nb = npart_ - na;
//			s += perfectIndexPartial(ketA)*comb_(nsite_,nb);
//			s += perfectIndexPartial(ketB);
//			if (s>=data_.size())
//				throw std::runtime_error("BasisOneSpinFeAs::PerfectIndex>=data_.size()\n");
//			return s;
//		}

		size_t getN(WordType ket,size_t site,size_t orb) const
		{
			std::vector<WordType> kets(orbitals_,0);
			uncollateKet(kets,ket);

			WordType res = (kets[orb] & bitmask_[site]);
			return (res>0) ? 1 : 0;
		}

		size_t getN(size_t i,size_t orb) const
		{
			std::vector<WordType> kets(orbitals_,0);
			uncollateKet(kets,data_[i]);
			return PsimagLite::BitManip::count(kets[orb]);
		}

		size_t getN(size_t i) const
		{
			size_t c = 0;
			for (size_t orb=0;orb<orbitals_;orb++)
				c += getN(i,orb);
			return c;
		}

		bool getBra(WordType& bra,const WordType& myword,size_t what,size_t site,size_t orb) const
		{
			std::vector<WordType> kets(orbitals_,0);
			uncollateKet(kets,myword);

			WordType braA = kets[orb];
			if (!getBra(braA,kets[orb],what,site)) return false;

			kets[orb] = braA;
			bra = getCollatedKet(kets);
			return true;
		}

		static const WordType& bitmask(size_t i)
		{
			return bitmask_[i];
		}

		int doSign(size_t i,size_t site,size_t orb) const
		{
			std::vector<WordType> kets(orbitals_,0);
			uncollateKet(kets,data_[i]);

			size_t c = 0;
			for (size_t orb1=0;orb1<orb;orb1++) {
				c += PsimagLite::BitManip::count(kets[orb1]);
			}

			int ret = (c&1) ? FERMION_SIGN : 1;
			return ret * doSign(kets[orb],site);
		}

		int doSign(
				WordType ket,
				size_t i,
				size_t orb1,
				size_t j,
				size_t orb2) const
		{
			if (i > j) {
				std::cerr<<"FATAL: At doSign\n";
				std::cerr<<"INFO: i="<<i<<" j="<<j<<std::endl;
				std::cerr<<"AT: "<<__FILE__<<" : "<<__LINE__<<std::endl;
				throw std::runtime_error("FeBasedSc::doSign(...)\n");
			}
			size_t x0 = (i+1)*orbitals_; // i+1 cannot be the last site, 'cause i<j
			size_t x1 = j*orbitals_;

			size_t sum = getNbyKet(ket,x0,x1);

			// at site i we need to be carefull
			x0 = i*orbitals_+orb1;
			x1 = (i+1)*orbitals_;
			sum += getNbyKet(ket,x0,x1);

			// same at site j
			x0 = j*orbitals_;
			x1 = j*orbitals_+orb2;
			sum += getNbyKet(ket,x0,x1);

			return (sum & 1) ? FERMION_SIGN : 1;
		}

		size_t getNbyKet(size_t ket) const
		{
			size_t sum = 0;
			WordType ketCopy = ket;
			while(ketCopy) {
				if (ketCopy & 1) sum++;
				ketCopy <<= 1;
			}
			return sum;
		}

		size_t isThereAnElectronAt(size_t ket,size_t site,size_t orb) const
		{
			size_t x = site*orbitals_ + orb;
			return (ket & bitmask_[x]) ? 1 : 0;
		}

		size_t electrons() const { return npart_; }

		int newPart(size_t type,size_t orb) const
		{
			int newPart1=npart_;

			int c = (type&1) ? -1 : 1;
			newPart1 += c;

			if (newPart1<0) return -1;

			if (size_t(newPart1)>orbitals_*nsite_) return -1;

			return newPart1;
		}

		int doSignGf(WordType a,size_t ind,size_t orb) const
		{
			std::vector<WordType> kets(orbitals_,0);
			uncollateKet(kets,a);

			size_t c = 0;
			for (size_t orb1=0;orb1<orb;orb1++) {
				c += PsimagLite::BitManip::count(kets[orb1]);
			}
			int ret = (c&1) ? FERMION_SIGN : 1;

			return ret * doSignGf(kets[orb],ind);
		}

	private:

		int doSignGf(WordType b,size_t ind) const
		{
			if (ind==0) return 1;

			// ind>0 from now on
			size_t i = 0;
			size_t j = ind;
			WordType mask = b;
			mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
			int s=(PsimagLite::BitManip::count(mask) & 1) ? -1 : 1; // Parity of up between i and j
			// Is there a down at i?
			if (bitmask_[i] & b) s = -s;
			return s;
		}

		void fillPartialBasis(std::vector<WordType>& partialBasis,size_t npart)
		{
			/* compute size of basis */
			size_t hilbert=1;
			int n=nsite_;
			size_t m=1;
			for (;m<=npart;n--,m++)
				hilbert=hilbert*n/m;

			if (partialBasis.size()!=hilbert) {
				partialBasis.clear();
				partialBasis.resize(hilbert);
			}

			if (npart==0) {
				partialBasis[0]=0;
				return;
			}
			/* define basis states */
			WordType ket = (1ul<<npart)-1;
			for (size_t i=0;i<hilbert;i++) {
				partialBasis[i] = ket;
				n=m=0;
				for (;(ket&3)!=1;n++,ket>>=1) {
					m += ket&1;
				}
				ket = ((ket+1)<<n) ^ ((1<<m)-1);
			}
		}

		void collateBasis(size_t& counter,const std::vector<std::vector<WordType> >& basisA)
		{
			size_t total = 1;
			for (size_t orb=0;orb<orbitals_;orb++) {
				total *= basisA[orb].size();
			}

			for (size_t i=0;i<total;i++) {
				std::vector<WordType> kets(orbitals_);
				getKets(kets,basisA,i);
				WordType ket = getCollatedKet(kets);
				data_[counter++] = ket;
			}
		}

		void getKets(std::vector<WordType>& kets,const std::vector<std::vector<WordType> >& basisA,size_t ind) const
		{
			// ind = i0 + i1 * size0 + i2 * size0 * size1 + ...
			size_t tmp = ind;

			size_t sizes = 1;
			for (size_t orb=0;orb<orbitals_-1;orb++)
				sizes *= basisA[orb].size();

			for (size_t orb=1;orb<orbitals_;orb++) {
				size_t ix = size_t(tmp / sizes);
				tmp = ind % sizes;
				kets[orbitals_-orb] = basisA[orbitals_-orb][ix];
				sizes /= basisA[orbitals_-orb-1].size();
			}
			kets[0] = basisA[0][tmp];
		}

		void doCombinatorial()
		{
			/* look-up table for binomial coefficients */
			comb_.reset(orbitals_*nsite_+1,orbitals_*nsite_+1);

			for (size_t n=0;n<comb_.n_row();n++)
				for (size_t i=0;i<comb_.n_col();i++)
					comb_(n,i)=0;

			for (size_t n=0;n<comb_.n_row();n++) {
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
			bitmask_.resize(nsite_*orbitals_*orbitals_+1);
			bitmask_[0]=1ul;
			for (size_t i=1;i<bitmask_.size();i++)
				bitmask_[i] = bitmask_[i-1]<<1;
		}

		size_t perfectIndexPartial(WordType state) const
		{
			size_t n=0;
			for (size_t b=0,c=1;state>0;b++,state>>=1)
				if (state&1) n += comb_(b,c++);

			return n;
		}

		WordType getCollatedKet(const std::vector<WordType>& kets) const
		{
			size_t counter = 0;
			WordType ket = 0;

			size_t b = orAll(kets);

			std::vector<WordType> remA = kets;

			while(b) {
				for (size_t orb=0;orb<kets.size();orb++) {
					size_t bitA = (remA[orb] & 1);
					if (bitA) ket |=bitmask_[counter+orb];
					counter++;
					if (remA[orb]) remA[orb] >>= 1;
				}
			}
			return ket;
		}

		size_t orAll(const std::vector<WordType>& kets) const
		{
			size_t b = 0;
			for (size_t orb=0;orb<kets.size();orb++) b |= kets[orb];
			return b;
		}

		//! kets must be all zero here
		void uncollateKet(std::vector<WordType>& kets,WordType ket) const
		{
			size_t counter = 0;

			while(ket) {
				for (size_t orb=0;orb<kets.size();orb++) {
					size_t mask = (1<<orb);
					size_t bitA = (ket & mask);

					if (bitA) kets[orb] |= bitmask_[counter];
				}
				counter++;
				ket >>= orbitals_;
			}
		}

		bool getBra(WordType& bra, const WordType& ket,size_t what,size_t i) const
		{

			WordType si=(ket & bitmask_[i]);
			if (what==DESTRUCTOR) {
				if (si>0) {
					bra = (ket ^ bitmask_[i]);
				} else {
					return false; // cannot destroy, there's nothing
				}
			} else {
				if (si==0) {
					bra = (ket ^ bitmask_[i]);
				} else {
					return false; // cannot construct, there's already one
				}
			}
			return true;
		}

		int doSign(WordType a, size_t i) const
		{
			if (i==nsite_-1) return 1;

			a &= ((1 << (i+1)) - 1) ^ ((1 << nsite_) - 1);
			// Parity of single occupied between i and nsite-1
			int s=(PsimagLite::BitManip::count(a) & 1) ? FERMION_SIGN : 1;
			return s;
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

		size_t size_;
		size_t npart_;
		std::vector<WordType> data_;

	}; // class BasisOneSpinFeAs

	std::ostream& operator<<(std::ostream& os,const BasisOneSpinFeAs& b)
	{
		for (size_t i=0;i<b.size();i++)
			os<<i<<" "<<b[i]<<"\n";
		return os;
	}

	size_t BasisOneSpinFeAs::orbitals_=2;
	size_t BasisOneSpinFeAs::nsite_=0;
	PsimagLite::Matrix<size_t> BasisOneSpinFeAs::comb_;
	std::vector<BasisOneSpinFeAs::WordType> BasisOneSpinFeAs::bitmask_;

} // namespace LanczosPlusPlus
#endif

