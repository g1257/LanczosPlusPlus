
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

#ifndef BASIS_ONE_SPIN_H
#define BASIS_ONE_SPIN_H
#include "Matrix.h"
#include "BitManip.h"

namespace LanczosPlusPlus {
	
	class BasisOneSpin {
	public:
		
		typedef unsigned int long long WordType;
		enum {DESTRUCTOR,CONSTRUCTOR};
		static size_t const ORBITALS  = 2;
		static int const FERMION_SIGN  = -1;
		static size_t nsite_;
		static PsimagLite::Matrix<size_t> comb_;
		static std::vector<WordType> bitmask_; 
		
		BasisOneSpin(size_t nsite, size_t npart)
				: npart_(npart)
		{
			if (nsite_>0 && nsite!=nsite_)
				throw std::runtime_error("BasisOneSpin: All basis must have same number of sites\n");
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
			for (size_t na=0;na<=npart;na++) {
				size_t nb = npart - na;
				std::vector<WordType> basisA, basisB;
				fillPartialBasis(basisA,na);
				fillPartialBasis(basisB,nb);
				collateBasis(counter,basisA,basisB);
			}
		}

		size_t size() const { return size_; }

		const WordType& operator[](size_t i) const
		{
			return data_[i];
		}

		size_t perfectIndex(WordType ket) const
		{
		//	for (size_t i=0;i<data_.size();i++)
		//		if (data_[i]==ket) return i;
		//	throw std::runtime_error("perfectindex\n");

			WordType ketA=0,ketB=0;
			uncollateKet(ketA,ketB,ket);
			// p(ket) = \sum_{na'=0}^{na'<na} S_na' * S_nb'
			//			+ p_A(ket_A)*S_nb + p_B(ket_B)
			// where S_x = C^n_x
			size_t na = PsimagLite::BitManip::count(ketA);
			// note nb = PsimagLite::BitManip::count(ketB)
			// or nb  = npart -na
			size_t s = 0;
			for (size_t nap=0;nap<na;nap++) {
				size_t nbp = npart_ - nap;
				s += comb_(nsite_,nap) * comb_(nsite_,nbp);
			}
			size_t nb = npart_ - na;
			s += perfectIndexPartial(ketA)*comb_(nsite_,nb);
			s += perfectIndexPartial(ketB);
			if (s>=data_.size())
				throw std::runtime_error("BasisOneSpin::PerfectIndex>=data_.size()\n");
			return s;
		}

		size_t getN(size_t i,size_t orb) const
		{
			WordType ketA=0,ketB=0;
			uncollateKet(ketA,ketB,data_[i]);
			if (orb==0) return PsimagLite::BitManip::count(ketA);
			return PsimagLite::BitManip::count(ketB);
		}

		size_t getN(size_t i) const
		{
			size_t c = 0;
			for (size_t orb=0;orb<ORBITALS;orb++)
				c += getN(i,orb);
			return c;
		}

		size_t getBraIndex(size_t i,size_t what,size_t orb) const
		{
			WordType ketA=0,ketB=0;
			uncollateKet(ketA,ketB,data_[i]);
			WordType braA = ketA;
			WordType braB = ketB;
			if (orb==0) {
				if (!getBra(braA,ketA,what,i)) throw std::runtime_error("getBraIndex::orb==0 problem\n");
			} else {
				if (!getBra(braB,ketB,what,i)) throw std::runtime_error("getBraIndex::orb!=0 problem\n");
			}
			WordType bra = getCollatedKet(braA,braB);
			return perfectIndex(bra);
		}

		static const WordType& bitmask(size_t i)
		{
			return bitmask_[i];
		}

		int doSign(size_t i,size_t site,size_t orb) const
		{
			WordType ketA=0,ketB=0;
			uncollateKet(ketA,ketB,data_[i]);
			if (orb==0) {
				return doSign(ketA,site);
			}

			size_t c = PsimagLite::BitManip::count(ketA);
			int ret = (c&1) ? FERMION_SIGN : 1;
			return ret * doSign(ketB,site);
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
			size_t x0 = (i+1)*ORBITALS; // i+1 cannot be the last site, 'cause i<j
			size_t x1 = j*ORBITALS;

			size_t sum = getNbyKet(ket,x0,x1);

			// at site i we need to be carefull
			x0 = i*ORBITALS+orb1;
			x1 = (i+1)*ORBITALS;
			sum += getNbyKet(ket,x0,x1);

			// same at site j
			x0 = j*ORBITALS;
			x1 = j*ORBITALS+orb2;
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
			size_t x = site*ORBITALS + orb;
			return (ket & bitmask_[x]) ? 1 : 0;
		}

	private:

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

		void collateBasis(
				size_t& counter,
				const std::vector<WordType>& basisA,
				const std::vector<WordType>& basisB)
		{
			for (size_t i=0;i<basisA.size();i++) {
				for (size_t j=0;j<basisB.size();j++) {
					WordType ket = getCollatedKet(basisA[i],basisB[j]);
					data_[counter++] = ket;
				}
			}
		}

		void doCombinatorial()
		{
			/* look-up table for binomial coefficients */
			comb_.reset(nsite_+1,nsite_+1);

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
			bitmask_.resize(nsite_*4+1);
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

		WordType getCollatedKet(WordType ketA,WordType ketB) const
		{
			WordType remA = ketA;
			WordType remB = ketB;
			size_t counter = 0;
			WordType ket = 0;

			while(remA || remB) {
				size_t bitA = (remA & 1);
				size_t bitB = (remB & 1);
				if (bitA) ket |=bitmask_[counter];
				if (bitB)  ket |=bitmask_[counter+1];
				counter += 2;
				if (remA) remA >>= 1;
				if (remB) remB >>= 1;
			}
			return ket;
		}

		void uncollateKet(WordType& ketA,WordType& ketB,WordType ket) const
		{
			size_t counter = 0;
			ketA = ketB = 0;
			while(ket) {
				size_t bitA = (ket & 1);
				size_t bitB = (ket & 2);
				if (bitA) ketA |= bitmask_[counter];
				if (bitB) ketB |= bitmask_[counter];
				counter++;
				ket >>= 2;
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

	}; // class BasisOneSpin

	std::ostream& operator<<(std::ostream& os,const BasisOneSpin& b)
	{
		for (size_t i=0;i<b.size();i++)
			os<<i<<" "<<b[i]<<"\n";
		return os;
	}

	size_t BasisOneSpin::nsite_=0;
	PsimagLite::Matrix<size_t> BasisOneSpin::comb_;
	std::vector<BasisOneSpin::WordType> BasisOneSpin::bitmask_;

} // namespace LanczosPlusPlus
#endif

