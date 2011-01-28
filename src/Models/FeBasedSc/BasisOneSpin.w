\documentclass{report}
\usepackage[T1]{fontenc}
\usepackage{bera}

\usepackage[pdftex,usenames,dvipsnames]{color}
\usepackage{listings}
\definecolor{mycode}{rgb}{0.9,0.9,1}
\lstset{language=c++,tabsize=1,basicstyle=\scriptsize,backgroundcolor=\color{mycode}}

\usepackage{hyperref}
\usepackage{fancyhdr}
\pagestyle{fancy}
\usepackage{verbatim}
\begin{document}

\begin{comment}
@o BasisOneSpin.h -t
@{
/*
@i license.txt
*/
@}
\end{comment}

\section{Basis a FeAsBasedSc Model to use with Lanczos}

HEre is some boilerplate:

@o BasisOneSpin.h -t
@{
#ifndef BASIS_ONE_SPIN_H
#define BASIS_ONE_SPIN_H
#include "Matrix.h"

namespace LanczosPlusPlus {
	@<theClassHere@>
} // namespace LanczosPlusPlus
#endif

@}

And the class is:
@d theClassHere
@{
class BasisOneSpin {
public:
	@<publicTypesAndEnums@>
	@<constructor@>
	@<publicFunctions@>

private:
	@<privateFunctions@>
	@<privateData@>
}; // class BasisOneSpin
@<staticDefinitions@>
@}

@d privateData
@{
size_t size_;
size_t npart_;
std::vector<WordType> data_;
@}

@d publicTypesAndEnums
@{
typedef unsigned int long long WordType;
static size_t nsite_;
static PsimagLite::Matrix<size_t> comb_;
static std::vector<WordType> bitmask_; @}

All right, now the constructor:
@d constructor
@{
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
		data_[0]=0;
		return;
	}
	size_ = 0;
	for (size_t na=0;na<=npart) {
			nb = npart - na;
			size_ += comb_(n,na) * comb_(n,nb);
	}
	data_.resize(size_);

	// compute basis:
	for (size_t na=0;na<=npart) {
		nb = npart - na;
		fillPartialBasis(basisA,na);
		fillPartialBasis(basisB,nb);
		collateBasis(basisA,basisB);
	}
}
@}

@d publicFunctions
@{
@<size@>
@<operatorBracket@>
@<perfectIndex@>
@<perfectIndexPartial@>
@<bitmask@>
@}

@d size
@{
size_t size() const { return size_; }
@}

@d operatorBracket
@{
const WordType& operator[](size_t i) const
{
	return data_[i];
}
@}

@d perfectIndex
@{
size_t perfectIndex(WordType ket) const
{
	WordType ketA,ketB;
	uncollateKet(ketA,ketB,ket);
	// p(ket) = \sum_{na'=0}^{na'<na} S_na' * S_nb'
	//			+ p_A(ket_A)*S_nb + p_B(ket_B)
	// where S_x = C^n_x
	size_t na = PsimagLite::BitManip::count(ketA);
	// note nb = PsimagLite::BitManip::count(ketB)
	// or nb  = npart -na
	size_t s = 0;
	for (size_t nap=0;nap<na;nap++) {
		nbp = npart - nap;
		s += comb_(nsite_,nap) * comb_(nsite_,nbp);
	}
	size_t nb = npart - na;
	s += perfectIndexPartial(ketA)*comb_(nsite,nb);
	s += perfectIndexPartial(ketB);
	return s;
}

@}

@d perfectIndexPartial
@{
size_t perfectIndexPartial(WordType state) const
{
	size_t n=0;
	for (size_t b=0,c=1;state>0;b++,state>>=1)
		if (state&1) n += comb_(b,c++);

	return n;
} @}

@d bitmask
@{
static const WordType& bitmask(size_t i)
{
	return bitmask_[i];
} @}


@d privateFunctions
@{
@<fillPartialBasis@>
@<collateBasis@>
@<doCombinatorial@>
@<doBitMask@>
@} 

@d fillPartialBasis
@{
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
@}

@d collateBasis
@{
void collateBasis()
{
	for (size_t i=0;i<basisA.size();i++) {
		for (size_t j=0;j<basisB.size();j++) {
			WordType ket = getCollatedKet(basisA[i],basisB[j]);
			data_[counter++] = ket;
		}
	}
}
@}

@d getCollatedKet
@{
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
@}

@d uncollateKet
@{
void uncollateKet(WordType& ketA,WordType& ketB,WordType ket) const
{
	while(ket) {
		size_t bitA = (ket & 1);
		size_t bitB = (ket & 2);
		if (bitA) ketA |= bitmask_[counter];
		if (bitB) ketB |= bitmask_[counter];
		counter++;
		ket >>= 2;
	}
}
@}


@d doCombinatorial
@{
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
@}

@d doBitMask
@{
void doBitmask()
{
	bitmask_.resize(nsite_);
	bitmask_[0]=1ul;
	for (size_t i=1;i<nsite_;i++)
		bitmask_[i] = bitmask_[i-1]<<1;
}
@}

@d staticDefinitions
@{
size_t BasisOneSpin::nsite_=0;
PsimagLite::Matrix<size_t> BasisOneSpin::comb_;
std::vector<typename BasisOneSpin::WordType> BasisOneSpin::bitmask_;
@}
	
\end{document}

