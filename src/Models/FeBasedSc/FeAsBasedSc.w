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
BasisType basis1_,basis2_;
@}

@d publicTypesAndEnums
@{
typedef BasisFeBasedSc BasisType;
typedef typename BasisType::WordType WordType;
typedef RealType_ RealType;
typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;
typedef std::vector<RealType> VectorType;
enum {SPIN_UP,SPIN_DOWN};
enum {DESTRUCTOR,CONSTRUCTOR};
static size_t const ORBITALS  = 1;
static size_t const DEGREES_OF_FREEDOM = 2*ORBITALS;
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
//void getOperator(SparseMatrixType& matrix,size_t what,size_t i,size_t sector) const
//{
//	size_t hilbert1 = basis1_.size();
//	size_t hilbert2 = basis2_.size();
//
//	matrix.resize(hilbert1*hilbert2);
//
//	size_t nCounter = 0;
//	for (size_t ispace1=0;ispace1<hilbert1;ispace1++) {
//		WordType ket1 = basis1_[ispace1];
//		for (size_t ispace2=0;ispace2<hilbert2;ispace2++) {
//			matrix.setRow(ispace2+hilbert2*ispace1,nCounter);
//
//			WordType ket2 = basis2_[ispace2];
//			WordType bra1 = ket1;
//			WordType bra2 = ket2;
//			if (sector==SPIN_DOWN) {
//				// modify bra1
//				if (!getBra(bra1,ket1,what,i)) continue;
//			} else {
//				//modify bra2:
//				if (!getBra(bra2,ket2,what,i)) continue;
//			}
//			size_t temp = perfectIndex(basis1_,basis2_,bra1,bra2);
//			matrix.pushCol(temp);
//			RealType cTemp=doSign(ket1,ket2,i,sector); // check SIGN FIXME
//
//			matrix.pushValue(cTemp);
//			nCounter++;
//		}
//	}
//	matrix.setRow(hilbert1*hilbert2,nCounter);
//}
@}

@d privateFunctions
@{
@<getBra@>
@<hoppings@>
@<setupHamiltonianP@>
@<countNonZero@>
@<perfectIndex@>
@<doSign1@>
@<doSign2@>
@}

This functions computes $|bra\rangle=c_i|ket\rangle$, where $c_i$ is a creation or destruction
operator (depending on whether \verb=what= is \verb=DESTRUCTOR= or \verb=CREATOR=)
on site $i$:
@d getBra
@{
//bool getBra(WordType& bra, const WordType& ket,size_t what,size_t i,size_t orbital) const
//{
//	size_t ii = i*DEGREES_OF_FREEDOM+orbital;
//	WordType si=(ket & BasisType::bitmask(ii));
//	if (what==DESTRUCTOR) {
//		if (si>0) {
//			bra = (ket ^ BasisType::bitmask(ii));
//		} else {
//			return false; // cannot destroy, there's nothing
//		}
//	} else {
//		if (si==0) {
//			bra = (ket ^ BasisType::bitmask(ii));
//		} else {
//			return false; // cannot contruct, there's already one
//		}
//	}
//	return true;
//}
@}

@d setupHamiltonianP
@{
void setupHamiltonian(SparseMatrixType &matrix,const BasisType &basis) const
{
	// Calculate diagonal elements AND count non-zero matrix elements
	size_t hilbert1=basis.size();
	std::vector<FieldType> diag(hilbert);
	size_t nzero = countNonZero(diag,basis);
	
	size_t nsite = geometry_.numberOfSites();
	
	// Setup CRS matrix
	matrix.resize(hilbert,nzero);
	
	// Calculate off-diagonal elements AND store matrix
	size_t nCounter=0;
	matrix.setRow(0,0);
	for (size_t ispace=0;ispace<hilbert;ispace++) {
		WordType ket1 = basis(ispace,SPIN_UP);
		WordType ket2 = basis(ispace],SPIN_DOWN);
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
							size_t temp = perfectIndex(basis,bra1,ket2);
							matrix.setCol(nCounter,temp);
							cTemp=h*doSign(ket1,ket2,i,orb,j,orb2,SPIN_DOWN); // check SIGN FIXME

							matrix.setValues(nCounter,cTemp);
							nCounter++;
						}
						if (s2i+s2j==1) {
							WordType bra2= ket2 ^(BasisType::bitmask(ii)|BasisType::bitmask(jj));
							size_t temp = perfectIndex(basis,ket1,bra2);
							matrix.setCol(nCounter,temp);
							cTemp=h*doSign(ket1,ket2,i,orb,j,orb2,SPIN_UP); // Check SIGN FIXME
							matrix.setValues(nCounter,cTemp);
							nCounter++;
						}
					}
				}
			}
		}
		matrix.setRow(ispace,nCounter);
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
size_t countNonZero(std::vector<FieldType>& diag,const BasisType &basis) const
{
	size_t hilbert=basis.size();
	size_t nsite = geometry_.numberOfSites();

	// Calculate diagonal elements AND count non-zero matrix elements
	size_t nzero = 0;
	for (size_t ispace=0;ispace<hilbert;ispace++) {
		WordType ket1 = basis[ispace].up;
		WordType ket2 = basis[ispace].down;
		RealType s=0;
		for (size_t i=0;i<nsite;i++) {
			for (size_t orb=0;orb<ORBITALS;orb++) {
				WordType s1i= (ket1 & BasisType::bitmask(i*ORBITALS+orb));
				WordType s2i= (ket2 & BasisType::bitmask(i*ORBITALS+orb));
				if (s1i>0) s1i=1;
				if (s2i>0) s2i=1;
				for (size_t orb2=0;orb2<ORBITALS;orb2++) {
					// Hubbard term
					s += mp_.hubbardU[i] * getN(ket1,orb) * getN(ket2,orb2);
				}
				
				// Potential term
				s += mp_.potentialV[i]*(getN(ket1,orb) + getN(ket2,orb));
				
				// Hopping term (only count how many non-zero)
				for (size_t j=0;j<nsite;j++) {
					if (j<i) continue;
					for (size_t orb2=0;orb2<ORBITALS;orb2++) {
						RealType tmp = hoppings(i,orb1,j,orb2);
						if (tmp==0) continue;
						WordType s1j= (ket1 & BasisType::bitmask(j*ORBITALS+orb2));
						WordType s2j= (ket2 & BasisType::bitmask(j*ORBITALS+orb2));
						if (s1j>0) s1j=1;
						if (s2j>0) s2j=1;
						if (s1i+s1j==1) nzero++;
						if (s2i+s2j==1) nzero++;
				}
			}
			// cout<<"diag of ("<<ispace1<<","<<ispace2<<"): "<<s<<endl;
			diag[ispace]=s;
			nzero++;
		}
	}

	nzero++;
	return nzero;
}
@}

@d perfectIndex
@{
size_t perfectIndex(const BasisType& basis,WordType ket1,WordType ket2) const
{
	size_t i1 = basis.perfectIndex(ket1,SPIN_UP);
	size_t i2 = basis.perfectIndex(ket2,SPIN_UP);
	return i1 + i2*basis.size(SPIN_UP);
}
@}

@d doSign1
@{
int doSign(WordType ket,size_t i,size_t flavor,size_t j,size_t flavor2) const
{
	if (i > j) {
		std::cerr<<"FATAL: At doSign\n";
		std::cerr<<"INFO: i="<<i<<" j="<<j<<std::endl;
		std::cerr<<"AT: "<<__FILE__<<" : "<<__LINE__<<std::endl;
		throw std::runtime_error("FeBasedSc::doSign(...)\n");
	}

	throw std::runtime_error("Unimplemented doSign\n");
}
@}

@d doSign2
@{
int doSign(WordType ket, size_t i,size_t flavor) const
{
	throw std::runtime_error("Unimplemented doSign2\n");
}
@}
			
\end{document}

