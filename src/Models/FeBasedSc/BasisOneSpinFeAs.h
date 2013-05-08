
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
#include "ProgramGlobals.h"

namespace LanczosPlusPlus {
	
	class BasisOneSpinFeAs {

		typedef Partitions PartitionsType;

		enum {OPERATOR_NIL=ProgramGlobals::OPERATOR_NIL,
		      OPERATOR_C=ProgramGlobals::OPERATOR_C,
		      OPERATOR_SZ=ProgramGlobals::OPERATOR_SZ,
		      OPERATOR_CDAGGER=ProgramGlobals::OPERATOR_CDAGGER};

	public:
		
		typedef ProgramGlobals::WordType WordType;
		static size_t orbitals_;
		static int const FERMION_SIGN  = -1;
		static size_t nsite_;
		static PsimagLite::Matrix<size_t> comb_;
		static PsimagLite::Vector<WordType>::Type bitmask_; 
		
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
			PartitionsType partitions(npart,orbitals_);
			for (size_t i=0;i<partitions.size();i++) {
					const PsimagLite::Vector<size_t>::Type& na = partitions(i);
					size_t tmp = 1;
					for (size_t j=0;j<na.size();j++)
						tmp *= comb_(nsite_,na[j]);
					size_ += tmp;
			}
			data_.resize(size_);

			// compute basis:
			size_t counter = 0;

			for (size_t i=0;i<partitions.size();i++) {
				const PsimagLite::Vector<size_t>::Type& na = partitions(i);
				PsimagLite::Vector<PsimagLite::Vector<WordType>::Type>::Type basisA(orbitals_);
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

		size_t getN(WordType ket,size_t site,size_t orb) const
		{
			PsimagLite::Vector<WordType>::Type kets(orbitals_,0);
			uncollateKet(kets,ket);

			WordType res = (kets[orb] & bitmask_[site]);
			return (res>0) ? 1 : 0;
		}

		size_t getN(size_t i,size_t orb) const
		{
			PsimagLite::Vector<WordType>::Type kets(orbitals_,0);
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
			PsimagLite::Vector<WordType>::Type kets(orbitals_,0);
			uncollateKet(kets,myword);

			WordType braA = kets[orb];
			if (!getBraCorCdagger(braA,kets[orb],what,site)) return false;

			kets[orb] = braA;
			bra = getCollatedKet(kets);
			return true;
		}

		static const WordType& bitmask(size_t i)
		{
			return bitmask_[i];
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

		int newPartCorCdagger(size_t what,size_t orb) const
		{
			int newPart1=npart_;

			int c = (what==ProgramGlobals::OPERATOR_C) ? -1 : 1;
			newPart1 += c;

			if (newPart1<0) return -1;

			if (size_t(newPart1)>orbitals_*nsite_) return -1;

			return newPart1;
		}

		int hasNewPartsSplusOrSminus(int c,size_t orb) const
		{
			int newPart1=npart_;

			newPart1 += c;

			if (newPart1<0) return -1;

			if (size_t(newPart1)>orbitals_*nsite_) return -1;

			return newPart1;
		}

		int doSignGf(WordType a,size_t ind,size_t orb) const
		{
			size_t x0 = 0;
			size_t x1 = ind*orbitals_;
			size_t sum = getNbyKet(a,x0,x1);

			// at site ind we need to be carefull
			x0 = ind*orbitals_;
			x1 = ind*orbitals_+orb;
			sum += getNbyKet(a,x0,x1);

			return (sum & 1) ? FERMION_SIGN : 1;
		}

	private:

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

		void collateBasis(size_t& counter,const PsimagLite::Vector<PsimagLite::Vector<WordType>::Type >::Type& basisA)
		{
			size_t total = 1;
			for (size_t orb=0;orb<orbitals_;orb++) {
				total *= basisA[orb].size();
			}

			for (size_t i=0;i<total;i++) {
				PsimagLite::Vector<WordType>::Type kets(orbitals_);
				getKets(kets,basisA,i);
				WordType ket = getCollatedKet(kets);
				data_[counter++] = ket;
			}
		}

		void getKets(PsimagLite::Vector<WordType>::Type& kets,const PsimagLite::Vector<PsimagLite::Vector<WordType>::Type >::Type& basisA,size_t ind) const
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

		WordType getCollatedKet(const PsimagLite::Vector<WordType>::Type& kets) const
		{
			size_t counter = 0;
			WordType ket = 0;

			PsimagLite::Vector<WordType>::Type remA = kets;

			while( orAll(remA) ) {
				for (size_t orb=0;orb<kets.size();orb++) {
					size_t bitA = (remA[orb] & 1);
					if (bitA) ket |=bitmask_[counter];
					counter++;
					if (remA[orb]) remA[orb] >>= 1;
				}
			}
			return ket;
		}

		size_t orAll(const PsimagLite::Vector<WordType>::Type& kets) const
		{
			size_t b = 0;
			for (size_t orb=0;orb<kets.size();orb++) b |= kets[orb];
			return b;
		}

		//! kets must be all zero here
		void uncollateKet(PsimagLite::Vector<WordType>::Type& kets,WordType ket) const
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

		bool getBraCorCdagger(WordType& bra, const WordType& ket,size_t what,size_t i) const
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
			std::string str = ProgramGlobals::unknownOperator(what);
			throw std::runtime_error(str.c_str());
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
		PsimagLite::Vector<WordType>::Type data_;

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
	PsimagLite::Vector<BasisOneSpinFeAs::WordType>::Type BasisOneSpinFeAs::bitmask_;

} // namespace LanczosPlusPlus
#endif

