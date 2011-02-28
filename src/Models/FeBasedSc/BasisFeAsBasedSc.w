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
@}

@d privateData
@{
BasisType basis1_,basis2_;
@}

@d publicTypesAndEnums
@{
typedef BasisOneSpin BasisType;
typedef BasisType::WordType WordType;
enum {SPIN_UP,SPIN_DOWN};
static size_t const ORBITALS  = BasisType::ORBITALS;
static int const FERMION_SIGN = BasisType::FERMION_SIGN;
enum {DESTRUCTOR=BasisType::DESTRUCTOR,CONSTRUCTOR=BasisType::CONSTRUCTOR};
@}

All right, now the constructor:
@d constructor
@{
BasisFeAsBasedSc(size_t nsite, size_t nup,size_t ndown)
		: basis1_(nsite,nup),basis2_(nsite,ndown)
{
	std::cout<<"Basis1\n";
	std::cout<<basis1_;
	std::cout<<"Basis2\n";
	std::cout<<basis2_;
}
@}

@d publicFunctions
@{
@<bitmask@>
@<size@>
@<operatorBracket@>
@<perfectIndex@>
@<getN@>
@<getBraIndex@>
@<doSign@>
@}

@d bitmask
@{
static const WordType& bitmask(size_t i)
{
	return BasisType::bitmask(i);
}
@}

@d size
@{
size_t size() const { return basis1_.size()*basis2_.size(); }
@}

@d operatorBracket
@{
const WordType& operator()(size_t i,size_t spin) const
{
	size_t y = i/basis1_.size();
	size_t x = i%basis1_.size();
	return (spin==SPIN_UP) ? basis1_[x] : basis2_[y];
}
@}

@d perfectIndex
@{
size_t perfectIndex(WordType ket1,WordType ket2) const
{
	return basis1_.perfectIndex(ket1) + basis2_.perfectIndex(ket2)*basis1_.size();
}
@}

@d getN
@{
size_t getN(size_t i,size_t spin,size_t orb) const
{
	size_t y = i/basis1_.size();
	size_t x = i%basis1_.size();
	return (spin==SPIN_UP) ? basis1_.getN(x,orb) : basis2_.getN(y,orb);
}
@}

@d getBraIndex
@{
size_t getBraIndex(size_t i,size_t what,size_t sector) const
{
	size_t y = i/basis1_.size();
	size_t x = i%basis1_.size();
	size_t i1 = basis1_.perfectIndex(basis1_[x]);
	size_t i2 = basis2_.perfectIndex(basis2_[y]);
	size_t spin = sector/2;
	size_t orb = (sector & 1);
	if (spin==SPIN_UP) {
		i1 = basis1_.getBraIndex(x,what,orb);
	} else {
		i2 =  basis2_.getBraIndex(y,what,orb);
	}

	return i1 + i2*basis1_.size();
}
@}

@d doSign
@{
int doSign(size_t i,size_t site,size_t sector) const
{
	size_t y = i/basis1_.size();
	size_t x = i%basis1_.size();
	size_t spin = sector/2;
	size_t orb = (sector & 1);
	if (spin==SPIN_UP) {
		return basis1_.doSign(x,site,orb);
	}
	size_t c = basis1_.getN(x);
	int ret = 1;
	if (c&1) ret = FERMION_SIGN;
	return ret * basis2_.doSign(y,site,orb);
}
@}

@d privateFunctions
@{

@} 


\end{document}

