/*
*/

#ifndef HUBBARDRASHBASOC_H
#define HUBBARDRASHBASOC_H

#include "BasisRashbaSOC.h"
#include "BitManip.h"
#include "TypeToString.h"
#include "SparseRow.h"
#include "../HubbardOneOrbital/ParametersModelHubbard.h"
#include "ProgramGlobals.h"
#include "../../Engine/ModelBase.h"
#include "../HubbardOneOrbital/HubbardHelper.h"

namespace LanczosPlusPlus {

template<typename ComplexOrRealType,typename GeometryType,typename InputType>
class HubbardOneOrbitalRashbaSOC : public ModelBase<ComplexOrRealType,GeometryType,InputType> {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef ModelBase<ComplexOrRealType,GeometryType,InputType> BaseType;
	typedef ParametersModelHubbard<RealType,InputType> ParametersModelType;
	typedef BasisRashbaSOC<GeometryType> BasisType;
	typedef typename BasisType::PairIntType PairIntType;
	typedef typename BasisType::BaseType BasisBaseType;
	typedef typename BasisType::WordType WordType;
	typedef typename BaseType::VectorSizeType VectorSizeType;
	typedef typename BaseType::SparseMatrixType SparseMatrixType;
	typedef typename BaseType::VectorType VectorType;
	typedef typename BasisType::LabeledOperatorType LabeledOperatorType;
	typedef PsimagLite::SparseRow<SparseMatrixType> SparseRowType;
	typedef HubbardHelper<BaseType, BasisType, ParametersModelType> HubbardHelperType;

	static int const FERMION_SIGN = BasisType::FERMION_SIGN;

	enum class TermEnum {HOPPING, RASHBA_SOC};

	HubbardOneOrbitalRashbaSOC(SizeType ne,
	                           InputType& io,
	                           const GeometryType& geometry)
	    : mp_(io),
	      geometry_(geometry),
	      basis_(geometry, ne),
	      helper_(geometry, mp_)
	{}

	~HubbardOneOrbitalRashbaSOC()
	{
		BaseType::deleteGarbage(garbage_);
	}

	SizeType size() const { return basis_.size(); }

	SizeType orbitals(SizeType) const
	{
		return 1;
	}

	void setupHamiltonian(SparseMatrixType& matrix) const
	{
		setupHamiltonian(matrix, basis_);
	}

	//! Gf. related functions below:
	void setupHamiltonian(SparseMatrixType& matrix,
	                      const BasisBaseType& basis) const
	{
		return helper_.setupHamiltonian(matrix, basis);
	}

	void matrixVectorProduct(VectorType &x,VectorType const &y) const
	{
		matrixVectorProduct(x,y,basis_);
	}

	void matrixVectorProduct(VectorType &x,
	                         VectorType const &y,
	                         const BasisBaseType& basis) const
	{
		helper_.matrixVectorProduct(x, y, basis);
	}

	bool hasNewParts(std::pair<SizeType,SizeType>& newParts,
	                 const std::pair<SizeType,SizeType>& oldParts,
	                 const LabeledOperator& lOperator,
	                 SizeType spin,
	                 SizeType orb) const
	{
		throw std::runtime_error("hasNewParts unimplemented\n");
	}

	const GeometryType& geometry() const { return geometry_; }

	const BasisType& basis() const { return basis_; }

	PsimagLite::String name() const { return __FILE__; }

	BasisType* createBasis(SizeType nup, SizeType ndown) const
	{
		BasisType* ptr = new BasisType(geometry_, nup + ndown);
		garbage_.push_back(ptr);
		return ptr;
	}

	void print(std::ostream& os) const { os<<mp_; }

	void printOperators(std::ostream& os) const
	{
		throw std::runtime_error("printOperators unimplemented\n");
	}

private:

	const ParametersModelType mp_;
	const GeometryType& geometry_;
	BasisType basis_;
	HubbardHelperType helper_;
	mutable typename PsimagLite::Vector<BasisType*>::Type garbage_;
}; // class HubbardOneOrbitalRashbaSOC
} // namespace LanczosPlusPlus
#endif

