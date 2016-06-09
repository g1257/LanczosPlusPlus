/*
Copyright (c) 2014-2016, UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/

#ifndef LANCZOS_REDUCED_DM_H
#define LANCZOS_REDUCED_DM_H
#include "Vector.h"
#include "Matrix.h"

namespace LanczosPlusPlus {

template<typename ModelType>
class ReducedDensityMatrix {

	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::BasisBaseType::WordType WordType;
	typedef typename ModelType::ComplexOrRealType ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef std::pair<SizeType,SizeType> PairSizeType;

public:

	ReducedDensityMatrix(const ModelType& model,
	                     const VectorType& psi,
	                     SizeType split)
		: row_(pow(model.basis().hilbertOneSite(),split)),
	      nabits_(split),
	      nbbits_(model.geometry().numberOfSites() - split),
	      rdm_(row_,row_)
	{
		build(model, psi);
		w_ = rdm_;
		diag(w_,eigs_,'V');
	}

	void printAll(std::ostream& os) const
	{
		os<<"Reduced Density Matrix\n";
		os<<rdm_;
		os<<"Eigenvectors of Reduced Density Matrix\n";
		os<<w_;
		os<<"Eigenvalues of Reduced Density Matrix\n";
		os<<eigs_;
	}

private:

	void build(const ModelType& model, const VectorType& psi)
	{
		SizeType hilbert = model.basis().size();
		for (SizeType i = 0; i < hilbert; ++i) {
			PairSizeType alphaBeta = unpack(model,i);
			for (SizeType j = 0; j < hilbert; ++j) {
				PairSizeType alphaPBeta = unpack(model,j);
				if (alphaBeta.second != alphaPBeta.second) continue;
				rdm_(alphaBeta.first,alphaPBeta.first) += std::conj(psi[i])*psi[j];
			}
		}
	}

	PairSizeType unpack(const ModelType& model, SizeType ind) const
	{
		WordType a = model.basis()(ind,0);
		WordType b = a;
		WordType mask = (1<<nabits_) - 1;
		a &= mask;

		mask = (1<<nbbits_) - 1;
		mask <<= nabits_;
		b &= mask;
		b >>= nabits_;
		return PairSizeType(a,b);
	}

	SizeType row_;
	SizeType nabits_;
	SizeType nbbits_;
	MatrixType rdm_;
	MatrixType w_;
	VectorRealType eigs_;
}; // class ReducedDensityMatrix

} // namespace LanczosPlusPlus
#endif

