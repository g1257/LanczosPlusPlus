
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2014, UT-Battelle, LLC
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

#ifndef LANCZOS_MODEL_BASE_H
#define LANCZOS_MODEL_BASE_H
#include "CrsMatrix.h"
#include "BasisBase.h"
#include "Vector.h"
#include "RahulOperator.h"

namespace LanczosPlusPlus {

template<typename ComplexOrRealType_,typename GeometryType_,typename InputType_>

class ModelBase {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType_>::Type RealType;
	typedef GeometryType_ GeometryType;
	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef InputType_ InputType;
	typedef BasisBase<GeometryType> BasisBaseType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef RahulOperator<ComplexOrRealType> RahulOperatorType;
	typedef typename PsimagLite::Vector<RahulOperatorType>::Type VectorRahulOperatorType;

	virtual ~ModelBase() {}

	virtual SizeType size() const = 0;

	virtual SizeType orbitals(SizeType site) const = 0;

	virtual void setupHamiltonian(SparseMatrixType& matrix) const = 0;

	virtual bool hasNewParts(std::pair<SizeType,SizeType>& newParts,
	                         const std::pair<SizeType,SizeType>& oldParts,
	                         const LabeledOperator& what,
	                         SizeType spin,
	                         SizeType orb) const = 0;

	virtual const GeometryType& geometry() const = 0;

	virtual void setupHamiltonian(SparseMatrixType& matrix,
	                              const BasisBaseType& basis) const =0;

	virtual void matrixVectorProduct(VectorType&,const VectorType&) const
	{
		throw PsimagLite::RuntimeError
		        ("ModelBase::matrixVectorProduct(2) not impl. for this model\n");
	}

	virtual void matrixVectorProduct(VectorType&,
	                                 const VectorType&,
	                                 const BasisBaseType&) const
	{
		throw PsimagLite::RuntimeError
		        ("ModelBase::matrixVectorProduct(3) not impl. for this model\n");
	}

	virtual const BasisBaseType& basis() const = 0;

	virtual PsimagLite::String name() const  = 0;

	virtual BasisBaseType* createBasis(SizeType nup, SizeType ndown) const = 0;

	virtual void print(std::ostream& os) const = 0;

	virtual void rahulMethod(VectorType&,
	                         const VectorRahulOperatorType&,
	                         const VectorSizeType&,
	                         const VectorType&) const
	{
		throw PsimagLite::RuntimeError("ModelBase::rahulMethod() not impl. for this model\n");
	}

	virtual void printOperators(std::ostream&) const
	{
		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) +  "\n";
		str += PsimagLite::String("Function printOperators unimplemented\n");
		std::cerr<<str.c_str();
	}

protected:

	template<typename SomeVectorType>
	static typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,
	void>::Type deleteGarbage(SomeVectorType& garbage)
	{
		for (SizeType i = 0; i < garbage.size(); ++i) {
			delete garbage[i];
			garbage[i] = 0;
		}
	}

}; // class ModelBase

template<typename RealType,typename GeometryType,typename InputType>
std::ostream& operator<<(std::ostream& os,
                         const ModelBase<RealType,GeometryType,InputType>& model)
{
	model.print(os);
	return os;
}

} // namespace LanczosPlusPlus

#endif

