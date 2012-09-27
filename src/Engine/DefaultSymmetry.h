
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
		typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;
		typedef std::vector<RealType> VectorType;

	public:

		DefaultSymmetry(const BasisType& basis,const GeometryType& geometry)
		{}

		void transformMatrix(std::vector<SparseMatrixType>& matrix1,const SparseMatrixType& matrix) const
		{
			throw std::runtime_error("DefaultSymmetry: cannot call transformMatrix\n");
		}

		void transformGs(VectorType& gs,size_t offset)
		{
		}

		size_t sectors() const { return 1; }

		std::string name() const { return "default"; }

	}; // class DefaultSymmetry
} // namespace Dmrg

#endif  // DEFAULT_SYMM_H
