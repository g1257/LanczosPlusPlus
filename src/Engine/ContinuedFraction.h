
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

#ifndef CONTINUED_FRACTION_H 
#define CONTINUED_FRACTION_H
#include <iostream>
#include "ProgressIndicator.h"
#include "BLAS.h"
#include "LanczosSolver.h"
#include "Random48.h"
#include "ProgramGlobals.h"

namespace LanczosPlusPlus {
	template<
    	typename ModelType_,
	 	typename ConcurrencyType_>
	class ContinuedFraction  {
	public:
		
		typedef ModelType_ ModelType;
		typedef ConcurrencyType_ ConcurrencyType;
		typedef typename ModelType::RealType RealType;
		typedef typename std::complex<RealType> ComplexType;
		typedef typename ModelType::VectorType VectorType;
		typedef typename ModelType::SparseMatrixType SparseMatrixType;
		typedef typename VectorType::value_type FieldType;
		typedef typename ModelType::BasisType BasisType;
		typedef PsimagLite::Random48<RealType> RandomType;
		typedef PsimagLite::LanczosSolver<RealType,SparseMatrixType,VectorType,RandomType,ProgramGlobals> LanczosSolverType;
		typedef PsimagLite::Matrix<FieldType> MatrixType;
		typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;

		static const size_t parallelRank_ = 0; // ContF needs to support concurrency FIXME
		enum {PLUS,MINUS};
		
		
		ContinuedFraction(const ModelType& model)
			: model_(model),progress_("ContinuedFraction",0)
		{
			// printHeader();
			// task 1: Compute Hamiltonian and
			// task 2: Compute ground state |phi>
			computeGroundState();
			
		} 
		

		RealType gsEnergy() const
		{
			return gsEnergy_;
		} 

		void getGreenFunction(TridiagonalMatrixType& ab,RealType& norma,
				size_t i,size_t j,size_t plusOrMinus) const
		{
			// task 3: compute |initVector> =\sum_x c_x|phi>, where
			// c_x are some operator
			VectorType initVector;
			computeInitVector(initVector,i,j,plusOrMinus);

			// task 4: tridiag H starting with |initVector>
			triDiagonalize(ab,initVector);

			norma = initVector*initVector;
		}
		

		ComplexType continuedFraction(ComplexType z,const TridiagonalMatrixType& ab) const
		{
			static MatrixType T;
			static bool firstcall = true;
			static std::vector<RealType> eigs(T.n_row());
			if (firstcall) {
				ab.buildDenseMatrix(T);
				diag(T,eigs,'V');
				firstcall = false;
			}
			ComplexType sum(0.0);
			for (size_t i=0;i<T.n_row();i++) {
				sum += T(i,0)*T(i,0)/(z-eigs[i]);
			}
			return sum;
		}

		
		
	
	private:
		

		void computeGroundState()
		{
			model_.setupHamiltonian(hamiltonian_);
			MatrixType fm;
			crsMatrixToFullMatrix(fm,hamiltonian_);
			bool verbose = true;
			//printNonZero(fm,std::cerr);
			if (!isHermitian(fm,verbose)) {
				//std::cerr<<fm;
				throw std::runtime_error("Hamiltonian non Hermitian\n");
			}
			//std::cerr<<hamiltonian_;
			std::cerr<<"Done setting up Hamiltonian\n";

			//fullDiag(fm);

			RealType eps= 0.01*ProgramGlobals::LanczosTolerance;
			size_t iter= ProgramGlobals::LanczosSteps;
			size_t parallelRank = 0;

			LanczosSolverType lanczosSolver(hamiltonian_,iter,eps,parallelRank);
			gsVector_.resize(hamiltonian_.rank());
			lanczosSolver.computeGroundState(gsEnergy_,gsVector_);
		} 

		void computeInitVector(VectorType& initVector,size_t i,size_t j,size_t plusOrMinus) const
		{
			initVector.resize(model_.size());
			VectorType tmpVector(initVector.size());
			size_t spin = ModelType::SPIN_UP;
			size_t destruction = ModelType::DESTRUCTOR;
			SparseMatrixType ci;
			model_.getOperator(ci,destruction,i,spin);
			SparseMatrixType cj;
			model_.getOperator(cj,destruction,j,spin);
			ci.matrixVectorProduct(tmpVector,gsVector_);
			cj.matrixVectorProduct(initVector,gsVector_);
			if (plusOrMinus == PLUS) initVector += tmpVector;
			else initVector -= tmpVector;
		} 

		void triDiagonalize(TridiagonalMatrixType& ab,const VectorType& initVector) const
		{
			// tridiagonalize starting with tmpVector = c^\dagger_i|gsVector>
			MatrixType V;

			RealType eps= 0.01*ProgramGlobals::LanczosTolerance;
			size_t iter= ProgramGlobals::LanczosSteps;
			size_t parallelRank = 0;

			LanczosSolverType lanczosSolver(hamiltonian_,iter,eps,parallelRank);

			lanczosSolver.tridiagonalDecomposition(initVector,ab,V);

		}
		

		void fullDiag(MatrixType& fm) const
		{
			std::vector<RealType> e(fm.n_row());
			diag(fm,e,'N');
			for (size_t i=0;i<e.size();i++)
				std::cout<<e[i]<<"\n";
		}
		
		
		
		const ModelType& model_;
		PsimagLite::ProgressIndicator progress_;
		SparseMatrixType hamiltonian_;
		RealType gsEnergy_;
		VectorType gsVector_; 
	}; // class ContinuedFraction
} // namespace Dmrg

#endif 
