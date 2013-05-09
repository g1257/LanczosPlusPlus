
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

#ifndef BASIS_ONE_SPIN_IMMM_H
#define BASIS_ONE_SPIN_IMMM_H
#include "Matrix.h"
#include "BitManip.h"
#include <cassert>
#include "ProgramGlobals.h"

namespace LanczosPlusPlus {

	class BasisOneSpinImmm {

	public:
		
		typedef ProgramGlobals::WordType WordType;

		enum {OPERATOR_NIL=ProgramGlobals::OPERATOR_NIL,
		      OPERATOR_C=ProgramGlobals::OPERATOR_C,
		      OPERATOR_SZ=ProgramGlobals::OPERATOR_SZ,
		      OPERATOR_CDAGGER=ProgramGlobals::OPERATOR_CDAGGER};

		static int const FERMION_SIGN  = -1;
		static size_t nsite_;
		static PsimagLite::Matrix<size_t> comb_;
		static PsimagLite::Vector<WordType>::Type bitmask_; 

		BasisOneSpinImmm(const PsimagLite::Vector<size_t>::Type& orbsPerSite, size_t npart)
		: orbsPerSite_(orbsPerSite),npart_(npart)
		{
			if (nsite_>0 && orbsPerSite.size()!=nsite_) {
				PsimagLite::String s =  
					"BasisOneSpinImmm: All basis must have same number of sites\n";
				throw std::runtime_error(s.c_str());
			}

			nsite_ = orbsPerSite.size();
			doCombinatorial();
			doBitmask();

			/* compute size of basis */
			if (npart==0) {
				data_.resize(1);
				data_[0]=0;
				return;
			}
			size_t levels = 0;
			for (size_t i=0;i<orbsPerSite_.size();i++)
				levels += orbsPerSite_[i];
			size_t tmp = comb_(levels,npart);
			data_.resize(tmp);
			
			levels = orbsPerSite_.size()*orbs();			
			tmp = comb_(levels,npart);
			reordering_.resize(tmp);

			// compute basis:
			size_t counter = 0;
			size_t counter2 = 0;
			for (size_t na=0;na<=npart;na++) {
				size_t nb = npart - na;
				PsimagLite::Vector<WordType>::Type basisA, basisB;
				fillPartialBasis(basisA,na);
				fillPartialBasis(basisB,nb);
				collateBasis(counter,counter2,basisA,basisB);
			}
// 			std::cerr<<" in ctor NPART="<<npart_<<"\n";
// 			print(std::cout);
		}

		size_t size() const { return data_.size(); }

		const WordType& operator[](size_t i) const
		{
			return data_[i];
		}
		
		void print(std::ostream& os) const
		{
			std::cerr<<"--------------npart="<<npart_<<"\n";
			for (size_t i=0;i<data_.size();i++)
				std::cerr<<data_[i]<<"\n";
			std::cerr<<"--------------\n";
		}

		size_t perfectIndex(WordType ket) const
		{
			for (size_t i=0;i<data_.size();i++)
				if (data_[i]==ket) return i;
			print(std::cout);
			assert(false);
			return 0; 
// 			WordType ketA=0,ketB=0;
// 			uncollateKet(ketA,ketB,ket);
// 			// p(ket) = \sum_{na'=0}^{na'<na} S_na' * S_nb'
// 			//			+ p_A(ket_A)*S_nb + p_B(ket_B)
// 			// where S_x = C^n_x
// 			size_t na = PsimagLite::BitManip::count(ketA);
// 			// note nb = PsimagLite::BitManip::count(ketB)
// 			// or nb  = npart -na
// 			size_t s = 0;
// 			for (size_t nap=0;nap<na;nap++) {
// 				size_t nbp = npart_ - nap;
// 				s += comb_(nsite_,nap) * comb_(nsite_,nbp);
// 			}
// 			size_t nb = npart_ - na;
// 			s += perfectIndexPartial(ketA)*comb_(nsite_,nb);
// 			s += perfectIndexPartial(ketB);
// 			assert(s<reordering_.size());
// 			assert(reordering_[s]<data_.size());
// 			return reordering_[s];
		}

		size_t getN(WordType ket,size_t site,size_t orb) const
		{
			WordType ketA=0,ketB=0;
			uncollateKet(ketA,ketB,ket);
			if (orb==0) {
				   WordType res = (ketA & bitmask_[site]);
				   return (res>0) ? 1 : 0;
			}
			WordType res2 = ketB & bitmask_[site];
			return (res2>0) ? 1 : 0;
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
// 			size_t c = 0;
// 			for (size_t orb=0;orb<orbsPerSite_[i];orb++)
// 				c += getN(i,orb);
// 			return c;
			throw std::runtime_error("getN\n");
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
				throw std::runtime_error("BasisOneSpinImmm::doSign(...)\n");
			}
			size_t x0 = (i+1)*orbs(); // i+1 cannot be the last site, 'cause i<j
			size_t x1 = j*orbs();

			size_t sum = getNbyKet(ket,x0,x1);

			// at site i we need to be carefull
			x0 = i*orbs()+orb1;
			x1 = (i+1)*orbs();
			sum += getNbyKet(ket,x0,x1);

			// same at site j
			x0 = j*orbs();
			x1 = j*orbs()+orb2;
			sum += getNbyKet(ket,x0,x1);

			return (sum & 1) ? FERMION_SIGN : 1;
		}
		
		int doSignGf(WordType a,size_t ind,size_t orb) const
		{
			WordType ketA=0,ketB=0;

			uncollateKet(ketA,ketB,a);

			if (orb==0) return doSignGf(ketA,ind);
			int s=(PsimagLite::BitManip::count(ketA) & 1) ? -1 : 1; // Parity of a

			return s * doSignGf(ketB,ind);
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
			size_t x = site*orbs() + orb;
			return (ket & bitmask_[x]) ? 1 : 0;
		}

		size_t electrons() const { return npart_; }

		bool getBra(WordType& bra,const WordType& myword,size_t what,size_t site,size_t orb) const
		{
			WordType ketA=0,ketB=0;
			uncollateKet(ketA,ketB,myword);
			WordType braA = ketA;
			WordType braB = ketB;

			if (orb==0) {
				if (!getBra(braA,ketA,what,site)) return false;
			} else {
				if (!getBra(braB,ketB,what,site)) return false;
			}
			bra = getCollatedKet(braA,braB);
			return true;
		}

		int newPartCorCdagger(size_t what,size_t orb) const
		{
			int newPart1=npart_;

			int c = (what==ProgramGlobals::OPERATOR_CDAGGER) ? 1 : -1;
			newPart1 += c;

			if (newPart1<0) return -1;

			if (size_t(newPart1)>maxElectrons()) return -1;

			return newPart1;
		}

	private:

		bool getBra(WordType& bra, const WordType& ket,size_t what,size_t i) const
		{

			WordType si=(ket & bitmask_[i]);
			if (what==OPERATOR_C) {
				if (si>0) {
					bra = (ket ^ bitmask_[i]);
					return true;
				} else {
					return false; // cannot destroy, there's nothing
				}
			} else if (what==OPERATOR_CDAGGER) {
				if (si==0) {
					bra = (ket ^ bitmask_[i]);
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

		size_t maxElectrons() const
		{
			size_t sum = 0;
			for (size_t i=0;i<orbsPerSite_.size();i++)
				sum += orbsPerSite_[i];
			return sum;
		}

		size_t orbs() const { return orbsPerSite_[0]; }

		void fillPartialBasis(PsimagLite::Vector<WordType>::Type& partialBasis,size_t npart)
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

		void collateBasis(size_t& counter,
		                  size_t& counter2,
		                  const PsimagLite::Vector<WordType>::Type& basisA,
		                  const PsimagLite::Vector<WordType>::Type& basisB)
		{
			for (size_t i=0;i<basisA.size();i++) {
				for (size_t j=0;j<basisB.size();j++) {
					reordering_[counter2++]=counter;
					if (isForbiddenSite(basisB[j])) continue;
					WordType ket = getCollatedKet(basisA[i],basisB[j]);
					assert(counter<data_.size());
					data_[counter++] = ket;
				}
			}
		}

		bool isForbiddenSite(const WordType& ket) const
		{
			for (size_t i=0;i<orbsPerSite_.size();i++) {
				if (orbsPerSite_[i]>1) continue;
				WordType mask = (1<<i);
				if (mask & ket) return true;
			}
			return false;
		}

		void doCombinatorial()
		{
			/* look-up table for binomial coefficients */
			comb_.reset(orbs()*nsite_+1,orbs()*nsite_+1);

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

		const PsimagLite::Vector<size_t>::Type& orbsPerSite_;
		size_t npart_;
		PsimagLite::Vector<WordType>::Type data_;
		PsimagLite::Vector<WordType>::Type reordering_;

	}; // class BasisOneSpinImmm

	std::ostream& operator<<(std::ostream& os,const BasisOneSpinImmm& b)
	{
		for (size_t i=0;i<b.size();i++)
			os<<i<<" "<<b[i]<<"\n";
		return os;
	}

	size_t BasisOneSpinImmm::nsite_=0;
	PsimagLite::Matrix<size_t> BasisOneSpinImmm::comb_;
	PsimagLite::Vector<BasisOneSpinImmm::WordType>::Type BasisOneSpinImmm::bitmask_;

} // namespace LanczosPlusPlus
#endif

