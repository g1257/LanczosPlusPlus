
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

#ifndef ENGINE_H_
#define ENGINE_H_
#include <iostream>
#include "ProgressIndicator.h"
#include "BLAS.h"
#include "LanczosSolver.h"
#include "Random48.h"
#include "ProgramGlobals.h"

namespace LanczosPlusPlus {
	template<typename ModelType_,typename InternalProductType,typename ConcurrencyType_>
	class Engine  {
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
		typedef PsimagLite::LanczosSolver<RealType,InternalProductType,
				VectorType,RandomType,ProgramGlobals> LanczosSolverType;
		typedef PsimagLite::Matrix<FieldType> MatrixType;
		typedef typename LanczosSolverType::TridiagonalMatrixType
				TridiagonalMatrixType;

		// ContF needs to support concurrency FIXME
		static const size_t parallelRank_ = 0;
		static const size_t CHECK_HERMICITY = 1;

		enum {PLUS,MINUS};
		
		Engine(const ModelType& model,size_t numberOfSites)
			: model_(model),numberOfSites_(numberOfSites),
			  progress_("ContinuedFraction",0)
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

		//! Calc Green function G(isite,jsite)  (still diagonal in spin)
		template<typename ContinuedFractionCollectionType>
		void greenFunction(
				ContinuedFractionCollectionType& cfCollection,
				int isite,
				int jsite,
				int spin) const
		{
			typedef typename ContinuedFractionCollectionType::
					ContinuedFractionType ContinuedFractionType;
			typedef typename ModelType::BasisType BasisType;

			for (size_t type=0;type<4;type++) {
				if (isite==jsite && type>1) continue;
				//if (type&1) continue;
				std::pair<size_t,size_t> newParts(0,0);
				if (!model_.hasNewParts(newParts,type,spin)) continue;
				// Create new bases
				BasisType basis1New(numberOfSites_,newParts.first);
				BasisType basis2New(numberOfSites_,newParts.second);

				std::vector<RealType> modifVector;
				model_.getModifiedState(modifVector,gsVector_,basis1New,basis2New,
						type,isite,jsite,spin);

				InternalProductType matrix(model_,basis1New,basis2New);
				ContinuedFractionType cf;

				calcGf(cf,modifVector,matrix,type,spin);
				cfCollection.push(cf);
			}
		}

	private:

		void computeGroundState()
		{
			InternalProductType hamiltonian(model_);
			//if (CHECK_HERMICITY) checkHermicity(h);

			RealType eps= ProgramGlobals::LanczosTolerance;
			size_t iter= ProgramGlobals::LanczosSteps;
			size_t parallelRank = 0;

			LanczosSolverType lanczosSolver(hamiltonian,iter,eps,parallelRank);
			gsVector_.resize(hamiltonian.rank());
			lanczosSolver.computeGroundState(gsEnergy_,gsVector_);
			std::cout<<"#GSNorm="<<(gsVector_*gsVector_)<<"\n";
		}

		template<typename ContinuedFractionType>
		void calcGf(ContinuedFractionType& cf,
				const std::vector<RealType>& modifVector,
				const InternalProductType& matrix,
				size_t type,
				size_t spin) const
		{
			typedef typename ContinuedFractionType::TridiagonalMatrixType
					TridiagonalMatrixType;

			RealType eps= ProgramGlobals::LanczosTolerance;
			size_t iter= ProgramGlobals::LanczosSteps;
			size_t parallelRank = 0;

			LanczosSolverType lanczosSolver(matrix,iter,eps,parallelRank);
			TridiagonalMatrixType ab;
			MatrixType V;
			lanczosSolver.tridiagonalDecomposition(modifVector,ab,V);
			RealType weight = modifVector*modifVector;
			//weight = 1.0/weight;
			int s = (type&1) ? -1 : 1;
			int s2 = (type>1) ? -1 : 1;
			//for (size_t i=0;i<ab.size();i++) ab.a(i) *= s;
			cf.set(ab,gsEnergy_,weight*s2,s);

		}
		
		//! For debugging purpose only:
		void fullDiag(MatrixType& fm) const
		{
			std::vector<RealType> e(fm.n_row());
			diag(fm,e,'N');
			for (size_t i=0;i<e.size();i++)
				std::cout<<e[i]<<"\n";
		}
		
		const ModelType& model_;
		size_t numberOfSites_;
		PsimagLite::ProgressIndicator progress_;
		RealType gsEnergy_;
		VectorType gsVector_; 
	}; // class ContinuedFraction
} // namespace Dmrg

#endif  // ENGINE_H_
