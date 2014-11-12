/*
Copyright (c) 2009-2014, UT-Battelle, LLC
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

#ifndef DEFAULT_SYMM_H
#define DEFAULT_SYMM_H
#include <iostream>
#include "ProgressIndicator.h"
#include "CrsMatrix.h"
#include "Vector.h"
#include "Matrix.h"

namespace LanczosPlusPlus {

template<typename GeometryType_,typename BasisType>
class DefaultSymmetry  {

	typedef typename GeometryType_::ComplexOrRealType ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef ProgramGlobals::WordType WordType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

public:

	typedef GeometryType_ GeometryType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

	DefaultSymmetry(const BasisType&,
	                const GeometryType&,
	                bool printMatrix)
	    : printMatrix_(printMatrix)
	{}

	template<typename SomeModelType>
	void init(const SomeModelType& model,const BasisType& basis)
	{
		model.setupHamiltonian(matrixStored_,basis);
		assert(isHermitian(matrixStored_));
		bool nrows = matrixStored_.row();
		if (printMatrix_) {
			if (nrows > 40)
				throw PsimagLite::RuntimeError("printMatrix too big\n");
			std::cout<<"#LanczosPlusPlus: DenseMatrix\n";
			std::cout<<matrixStored_.toDense();
			PsimagLite::Matrix<ComplexOrRealType> matrixCopy;
			VectorRealType eigs(matrixCopy.n_row());
			fullDiag(eigs,matrixCopy);
		}
	}

	void fullDiag(VectorRealType& eigs,MatrixType& fm) const
	{
		if (matrixStored_.row() > 1000)
			throw PsimagLite::RuntimeError("fullDiag too big\n");

		fm = matrixStored_.toDense();
		diag(fm,eigs,'V');

		if (!printMatrix_) return;

		for (SizeType i=0;i<eigs.size();i++)
			std::cout<<eigs[i]<<"\n";
		std::cout<<fm;
	}

	void transformMatrix(typename PsimagLite::Vector<SparseMatrixType>::Type& matrix1,
	                     const SparseMatrixType& matrix) const
	{
		throw std::runtime_error("DefaultSymmetry: cannot call transformMatrix\n");
	}

	void transformGs(VectorType&,SizeType)
	{
	}

	SizeType sectors() const { return 1; }

	void setPointer(SizeType) { }

	PsimagLite::String name() const { return "default"; }

	SizeType rank() const { return matrixStored_.row(); }

	template<typename SomeVectorType>
	void matrixVectorProduct(SomeVectorType &x, SomeVectorType const &y) const
	{
		return matrixStored_.matrixVectorProduct(x,y);
	}

private:

	SparseMatrixType matrixStored_;
	bool printMatrix_;
}; // class DefaultSymmetry
} // namespace Dmrg

#endif  // DEFAULT_SYMM_H

