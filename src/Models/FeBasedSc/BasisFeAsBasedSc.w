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
@o BasisFeAsBasedSc.h -t
@{
/*
@i license.txt
*/
@}
\end{comment}

\section{Basis a FeAsBasedSc Model to use with Lanczos}

HEre is some boilerplate:

@o BasisFeAsBasedSc.h -t
@{
#ifndef BASIS_FEASBASED_SC_H
#define BASIS_FEASBASED_SC_H
#include "BasisOneSpin.h"

namespace LanczosPlusPlus {
	@<theClassHere@>
} // namespace LanczosPlusPlus
#endif

@}

And the class is:
@d theClassHere
@{
class BasisFeAsBasedSc {
public:
	@<publicTypesAndEnums@>
	@<constructor@>
	@<publicFunctions@>

private:
	@<privateFunctions@>
	@<privateData@>
}; // class BasisFeAsBasedSc
@<staticDefinitions@>
@}

@d privateData
@{
BasisType basis1_,basis2_;
@}

@d publicTypesAndEnums
@{
typedef unsigned int long long WordType;
static size_t nsite_;
static PsimagLite::Matrix<size_t> comb_;
static std::vector<WordType> bitmask_;
@}

All right, now the constructor:
@d constructor
@{
BasisFeAsBasedSc(size_t nsite, size_t nup,size_t ndown)
		: basis1_(nsite,nup),basis2_(nsite,ndown)
{
}
@}

@d publicFunctions
@{
@<size@>
@<operatorBracket@>
@<perfectIndex@>
@<bitmask@>
@<bitctn@>
@}

@d size
@{
size_t size() const { return basis1_.size()*basis2_.size(); }
@}

@d operatorBracket
@{
//const WordType& operator[](size_t i) const
//{
//	return data_[i];
//}
@}

@d perfectIndex
@{
size_t perfectIndex(WordType state) const
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

@d bitctn
@{
static int bitcnt (WordType b)
{
#if (ULONG_MAX == 0xfffffffful)
	b = (0x55555555 & b) + (0x55555555 & (b >> 1));
	b = (0x33333333 & b) + (0x33333333 & (b >> 2));
	b = (0x0f0f0f0f & b) + (0x0f0f0f0f & (b >> 4));
	b = (0x00ff00ff & b) + (0x00ff00ff & (b >> 8));
	b = (0x0000ffff & b) + (0x0000ffff & (b >> 16));

	return (int) b;
#else
	b = (0x5555555555555555 & b) + (0x5555555555555555 & (b >> 1));
	b = (0x3333333333333333 & b) + (0x3333333333333333 & (b >> 2));
	b = (0x0f0f0f0f0f0f0f0f & b) + (0x0f0f0f0f0f0f0f0f & (b >> 4));
	b = (0x00ff00ff00ff00ff & b) + (0x00ff00ff00ff00ff & (b >> 8));
	b = (0x0000ffff0000ffff & b) + (0x0000ffff0000ffff & (b >> 16));
	b = (0x00000000ffffffff & b) + (0x00000000ffffffff & (b >> 32));

	return (int) b;
#endif
}
@}

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
WordType getCollatedKet()
{
	WordType remA = ketA;
	WordType remB = ketB;
	size_t counter = 0;
	WordType ket = 0;

	while(remA || remB) {
		size_t bitA = (remA & 1);
		size_t bitB = (remB & 1);
		if (bitA) {
			WordType a = bitmask_[counter];
			ket |= a;
		}
		if (bitB) {
			WordType b = bitmask_[counter+1];
			ket |= b;
		}
		counter += 2;
		remA >> 1;
		remB >> 1;
	}
	return ket;
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
size_t BasisFeAsBasedSc::nsite_=0;
PsimagLite::Matrix<size_t> BasisFeAsBasedSc::comb_;
std::vector<typename BasisFeAsBasedSc::WordType> BasisFeAsBasedSc::bitmask_;
@}
	
\end{document}

