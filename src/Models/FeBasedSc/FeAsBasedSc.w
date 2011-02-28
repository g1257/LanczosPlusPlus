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

%\title{The FeBasedSc Class}
%\author{G.A.}
%\maketitle

\begin{comment}
@o FeBasedSc.h -t
@{
/*
@i license.txt
*/
@}
\end{comment}

\section{A FeBasedSc model to use with Lanczos}

HEre is some boilerplate:

@o FeBasedSc.h -t
@{
#ifndef FEBASED_SC_H
#define FEBASED_SC_H

#include "CrsMatrix.h"
#include "BasisFeAsBasedSc.h"

namespace LanczosPlusPlus {
	@<theClassHere@>
} // namespace LanczosPlusPlus
#endif

@}

And the class is:
@d theClassHere
@{
template<typename RealType_,typename ParametersType,
	typename GeometryType>
class FeBasedSc {
	@<privateTypesAndEnums@>
public:
	@<publicTypesAndEnums@>
	@<constructor@>
	@<publicFunctions@>
private:
	@<privateFunctions@>
	@<privateData@>
}; // class FeBasedSc
@}

@d privateTypesAndEnums
@{
typedef PsimagLite::Matrix<RealType_> MatrixType;
@}

@d privateData
@{
const ParametersType& mp_;
const GeometryType& geometry_;
BasisType basis_;
@}

@d publicTypesAndEnums
@{
typedef BasisFeAsBasedSc BasisType;
typedef typename BasisType::WordType WordType;
typedef RealType_ RealType;
typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;
typedef std::vector<RealType> VectorType;
enum {SPIN_UP=BasisType::SPIN_UP,SPIN_DOWN=BasisType::SPIN_DOWN};
enum {DESTRUCTOR=BasisType::DESTRUCTOR,CONSTRUCTOR=BasisType::CONSTRUCTOR};
static size_t const ORBITALS  = BasisType::ORBITALS;
static size_t const DEGREES_OF_FREEDOM = 2*ORBITALS;
static int const FERMION_SIGN = BasisType::FERMION_SIGN;
@}

@d constructor
@{
FeBasedSc(size_t nup,size_t ndown,const ParametersType& mp,GeometryType& geometry)
	: mp_(mp),geometry_(geometry),
	  basis_(geometry.numberOfSites(),nup,ndown)
{
}
@}

@d publicFunctions
@{
@<size@>
@<setupHamiltonian@>
@<getOperator@>
@}

@d size
@{
	size_t size() const { return basis_.size(); }
@}

@d setupHamiltonian
@{
void setupHamiltonian(SparseMatrixType &matrix) const
{
	setupHamiltonian(matrix,basis_);
}
@}

@d getOperator
@{
void getOperator(SparseMatrixType& matrix,size_t what,size_t i,size_t flavor) const
{
	size_t hilbert = basis_.size();

	matrix.resize(hilbert);

	size_t nCounter = 0;
	for (size_t ispace=0;ispace<hilbert;ispace++) {
		matrix.setRow(ispace,nCounter);

		size_t temp = basis_.getBraIndex(ispace,what,flavor);
		matrix.pushCol(temp);
		RealType cTemp=basis_.doSign(ispace,i,flavor); // check SIGN FIXME

		matrix.pushValue(cTemp);
		nCounter++;
	}
	matrix.setRow(hilbert,nCounter);
}
@}

@d privateFunctions
@{
@<hoppings@>
@<setupHamiltonianP@>
@<countNonZero@>
@}

This functions computes $|bra\rangle=c_i|ket\rangle$, where $c_i$ is a creation or destruction
operator (depending on whether \verb=what= is \verb=DESTRUCTOR= or \verb=CREATOR=)
on site $i$:


@d setupHamiltonianP
@{
void setupHamiltonian(SparseMatrixType &matrix,const BasisType &basis) const
{
	// Calculate diagonal elements AND count non-zero matrix elements
	size_t hilbert=basis.size();
	std::vector<RealType> diag(hilbert);
	size_t nzero = countNonZero(diag,basis);
	
	size_t nsite = geometry_.numberOfSites();
	
	// Setup CRS matrix
	matrix.resize(hilbert,nzero);
	
	// Calculate off-diagonal elements AND store matrix
	size_t nCounter=0;
	matrix.setRow(0,0);
	for (size_t ispace=0;ispace<hilbert;ispace++) {
		matrix.setRow(ispace,nCounter);
		WordType ket1 = basis(ispace,SPIN_UP);
		WordType ket2 = basis(ispace,SPIN_DOWN);
		// Save diagonal
		matrix.setCol(nCounter,ispace);
		RealType cTemp=diag[ispace];
		matrix.setValues(nCounter,cTemp);
		nCounter++;
		for (size_t i=0;i<nsite;i++) {
			for (size_t orb=0;orb<ORBITALS;orb++) {
				size_t ii = i*ORBITALS+orb;
				WordType s1i=(ket1 & BasisType::bitmask(ii));
				if (s1i>0) s1i=1;
				WordType s2i=(ket2 & BasisType::bitmask(ii));
				if (s2i>0) s2i=1;

				// Hopping term 
				for (size_t j=0;j<nsite;j++) {
					if (j<i) continue;
					for (size_t orb2=0;orb2<ORBITALS;orb2++) {
						size_t jj = j*ORBITALS+orb2;
						RealType h = hoppings(i,orb,j,orb2);
						if (h==0) continue;
						WordType s1j= (ket1 & BasisType::bitmask(jj));
						if (s1j>0) s1j=1;
						WordType s2j= (ket2 & BasisType::bitmask(jj));
						if (s2j>0) s2j=1;

						if (s1i+s1j==1) {
							WordType bra1= ket1 ^(BasisType::bitmask(ii)|BasisType::bitmask(jj));
							size_t temp = basis.perfectIndex(bra1,ket2);
							matrix.setCol(nCounter,temp);
							int extraSign = (s1i==0) ? FERMION_SIGN : 1;
							cTemp=h*extraSign*basis_.doSign(
								ket1,ket2,i,orb,j,orb2,SPIN_UP); // check SIGN FIXME

							matrix.setValues(nCounter,cTemp);
							nCounter++;
						}
						if (s2i+s2j==1) {
							WordType bra2= ket2 ^(BasisType::bitmask(ii)|BasisType::bitmask(jj));
							size_t temp = basis.perfectIndex(ket1,bra2);
							matrix.setCol(nCounter,temp);
							int extraSign = (s2i==0) ? FERMION_SIGN : 1;
							cTemp=h*extraSign*basis_.doSign(
								ket1,ket2,i,orb,j,orb2,SPIN_DOWN); // Check SIGN FIXME
							matrix.setValues(nCounter,cTemp);
							nCounter++;
						}
					}
				}
			}
		}
	}
	matrix.setRow(hilbert,nCounter);
}
@}

@d hoppings
@{
RealType hoppings(size_t i,size_t orb1,size_t j,size_t orb2) const
{
	return geometry_(i,orb1,j,orb2,0);
}
@}
		
@d countNonZero
@{
size_t countNonZero(std::vector<RealType>& diag,const BasisType &basis) const
{
	size_t hilbert=basis.size();
	size_t nsite = geometry_.numberOfSites();

	// Calculate diagonal elements AND count non-zero matrix elements
	size_t nzero = 0;
	for (size_t ispace=0;ispace<hilbert;ispace++) {
		WordType ket1 = basis(ispace,SPIN_UP);
		WordType ket2 = basis(ispace,SPIN_DOWN);
		RealType s=0;
		for (size_t i=0;i<nsite;i++) {
			for (size_t orb=0;orb<ORBITALS;orb++) {
				WordType s1i= (ket1 & BasisType::bitmask(i*ORBITALS+orb));
				WordType s2i= (ket2 & BasisType::bitmask(i*ORBITALS+orb));
				if (s1i>0) s1i=1;
				if (s2i>0) s2i=1;
				for (size_t orb2=0;orb2<ORBITALS;orb2++) {
					// Hubbard term
					s += mp_.hubbardU[i] * basis.getN(ispace,SPIN_UP,orb) *
							basis.getN(ispace,SPIN_DOWN,orb2);
				}
				
				// Potential term
				s += mp_.potentialV[i]*(basis.getN(ispace,SPIN_UP,orb) *
						basis.getN(ispace,SPIN_DOWN,orb));
				
				// Hopping term (only count how many non-zero)
				for (size_t j=0;j<nsite;j++) {
					if (j<i) continue;
					for (size_t orb2=0;orb2<ORBITALS;orb2++) {
						RealType tmp = hoppings(i,orb,j,orb2);
						if (tmp==0) continue;
						WordType s1j= (ket1 & BasisType::bitmask(j*ORBITALS+orb2));
						WordType s2j= (ket2 & BasisType::bitmask(j*ORBITALS+orb2));
						if (s1j>0) s1j=1;
						if (s2j>0) s2j=1;
						if (s1i+s1j==1) nzero++;
						if (s2i+s2j==1) nzero++;
					}
				}
			}
		}
		// cout<<"diag of ("<<ispace1<<","<<ispace2<<"): "<<s<<endl;
		diag[ispace]=s;
		nzero++;
	}

	nzero++;
	return nzero;
}
@}
			
\end{document}

