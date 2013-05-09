
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

#ifndef DEFAULT_SYMM_H
#define DEFAULT_SYMM_H
#include <iostream>
#include "ProgressIndicator.h"
#include "CrsMatrix.h"
#include "Vector.h"

namespace LanczosPlusPlus {


	template<typename GeometryType,typename BasisType>
	class DefaultSymmetry  {

		typedef typename GeometryType::RealType RealType;
		typedef typename BasisType::WordType WordType;

	public:

		typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;
		typedef typename PsimagLite::Vector<RealType>::Type VectorType;

		DefaultSymmetry(const BasisType& basis,const GeometryType& geometry)
		{
		}

		template<typename SomeModelType>
		void init(const SomeModelType& model,const BasisType& basis)
		{
			model.setupHamiltonian(matrixStored_,basis);
//			std::cout<<matrixStored_;
		}

		void transformMatrix(typename PsimagLite::Vector<SparseMatrixType>::Type& matrix1,
		                     const SparseMatrixType& matrix) const
		{
			throw std::runtime_error("DefaultSymmetry: cannot call transformMatrix\n");
		}

		void transformGs(VectorType& gs,size_t offset)
		{
		}

		size_t sectors() const { return 1; }

		void setPointer(size_t p) { }

		PsimagLite::String name() const { return "default"; }

		size_t rank() const { return matrixStored_.row(); }

		template<typename SomeVectorType>
		void matrixVectorProduct(SomeVectorType &x, SomeVectorType const &y) const
		{
			return matrixStored_.matrixVectorProduct(x,y);
		}

	private:

		SparseMatrixType matrixStored_;

	}; // class DefaultSymmetry
} // namespace Dmrg

#endif  // DEFAULT_SYMM_H
