
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
#include "ParametersForSolver.h"
#include "ParametersEngine.h"

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
		typedef typename InternalProductType::ReflectionSymmetryType ReflectionSymmetryType;
		typedef PsimagLite::Random48<RealType> RandomType;
		typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
		typedef PsimagLite::LanczosSolver<ParametersForSolverType,InternalProductType,VectorType>
		                    LanczosSolverType;
		typedef PsimagLite::Matrix<FieldType> MatrixType;
		typedef typename LanczosSolverType::TridiagonalMatrixType
				TridiagonalMatrixType;

		// ContF needs to support concurrency FIXME
		static const size_t parallelRank_ = 0;
		static const size_t CHECK_HERMICITY = 1;

		enum {PLUS,MINUS};
		
		Engine(const ModelType& model,size_t numberOfSites,PsimagLite::IoSimple::In& io)
		: model_(model),
		  progress_("Engine",0),
		  params_(io)
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
		void spectralFunction(ContinuedFractionCollectionType& cfCollection,
							  size_t what,
							  int isite,
							  int jsite,
							  int spin,
							  const std::pair<size_t,size_t>& orbs) const
		{
			typedef typename ContinuedFractionCollectionType::ContinuedFractionType ContinuedFractionType;
			typedef typename ModelType::BasisType BasisType;
			const BasisType* basisNew = 0;

			for (size_t type=0;type<4;type++) {
				if (isite==jsite && type>1) continue;
				//if (type&1) continue;
				if (ProgramGlobals::needsNewBasis(what)) {
					std::pair<size_t,size_t> newParts(0,0);
					if (!model_.hasNewParts(newParts,what,type,spin,orbs)) continue;
					// Create new bases
					basisNew = new BasisType(model_.geometry(),newParts.first,newParts.second);
				} else {
					basisNew = &model_.basis();
				}
				std::vector<RealType> modifVector;
				model_.getModifiedState(modifVector,what,gsVector_,*basisNew,type,isite,jsite,spin);

				InternalProductType matrix(model_,*basisNew);
				ContinuedFractionType cf;

				calcSpectral(cf,modifVector,matrix,type,spin);
				cfCollection.push(cf);

				if (ProgramGlobals::needsNewBasis(what)) delete basisNew;
			}
		}

		void ciCj(PsimagLite::Matrix<RealType>& result,size_t spin,const std::pair<size_t,size_t>& orbs) const
		{
			size_t type = 0;
			std::pair<size_t,size_t> newParts(0,0);
			if (!model_.hasNewParts(newParts,ProgramGlobals::SPECTRAL_CC,type,spin,orbs)) return;

			BasisType basisNew(model_.geometry(),newParts.first,newParts.second);

			std::cerr<<"basisNew.size="<<basisNew.size()<<" ";
			std::cerr<<"newparts.first="<<newParts.first<<" ";
			std::cerr<<"newparts.second="<<newParts.second<<"\n";

			size_t isign = 1;

			size_t total =result.n_row();
			RealType sum = 0;
			for (size_t isite=0;isite<total;isite++) {
				std::vector<RealType> modifVector1(basisNew.size(),0);
				model_.accModifiedState(modifVector1,basisNew,gsVector_,BasisType::DESTRUCTOR,
							isite,spin,orbs.first,isign);
				for (size_t jsite=0;jsite<total;jsite++) {
					std::vector<RealType> modifVector2(basisNew.size(),0);
					model_.accModifiedState(modifVector2,basisNew,gsVector_,BasisType::DESTRUCTOR,
								jsite,spin,orbs.second,isign);
					result(isite,jsite) =  modifVector2*modifVector1;
				}
				sum += result(isite,isite);
			}
			std::cout<<"Total Electrons = "<<sum<<"\n";
		}

	private:

//		void computeGroundState()
//		{
//			InternalProductType hamiltonian(model_);
//			//if (CHECK_HERMICITY) checkHermicity(h);

//			RealType eps= ProgramGlobals::LanczosTolerance;
//			size_t iter= ProgramGlobals::LanczosSteps;

//			ParametersForSolverType params;
//			params.steps = iter;
//			params.tolerance = eps;
//			params.lotaMemory = lotaMemory_;
//			params.stepsForEnergyConvergence =ProgramGlobals::MaxLanczosSteps;

//			LanczosSolverType lanczosSolver(hamiltonian,params);

//			gsVector_.resize(hamiltonian.rank());
//			lanczosSolver.computeGroundState(gsEnergy_,gsVector_);
//			std::cout<<"#GSNorm="<<(gsVector_*gsVector_)<<"\n";
//		}

		void computeGroundState()
		{
			ReflectionSymmetryType* rs=0;
			if (params_.useReflectionSymmetry)
				rs = new ReflectionSymmetryType(model_.basis(),model_.geometry());
			InternalProductType hamiltonian(model_,rs);
			//if (CHECK_HERMICITY) checkHermicity(h);

			RealType eps= ProgramGlobals::LanczosTolerance;
			size_t iter= ProgramGlobals::LanczosSteps;

			ParametersForSolverType params;
			params.steps = iter;
			params.tolerance = eps;
			params.lotaMemory = params_.storeLanczosVectors;
			params.stepsForEnergyConvergence =ProgramGlobals::MaxLanczosSteps;

			LanczosSolverType lanczosSolver(hamiltonian,params);

			VectorType gsVector1(hamiltonian.rank());
			RealType gsEnergy1 = 0;
			lanczosSolver.computeGroundState(gsEnergy1,gsVector1);

			if (!params_.useReflectionSymmetry) {
				gsVector_=gsVector1;
				gsEnergy_=gsEnergy1;
				return;
			}

			hamiltonian.reflectionSector(1);
			VectorType gsVector2(hamiltonian.rank());
			RealType gsEnergy2 = 0;
			lanczosSolver.computeGroundState(gsEnergy2,gsVector2);

			gsEnergy_=rs->setGroundState(gsVector_,gsEnergy1,gsVector1,gsEnergy2,gsVector2);

			std::cout<<"#GSNorm="<<(gsVector_*gsVector_)<<"\n";
		}

		template<typename ContinuedFractionType>
		void calcSpectral(ContinuedFractionType& cf,
						  const std::vector<RealType>& modifVector,
						  const InternalProductType& matrix,
						  size_t type,
						  size_t spin) const
		{
			typedef typename ContinuedFractionType::TridiagonalMatrixType
			                                        TridiagonalMatrixType;
			

			RealType eps= ProgramGlobals::LanczosTolerance;
			size_t iter= ProgramGlobals::LanczosSteps;

			ParametersForSolverType params;
			params.steps = iter;
			params.tolerance = eps;
			params.lotaMemory = params_.storeLanczosVectors;
			params.stepsForEnergyConvergence =ProgramGlobals::MaxLanczosSteps;

			LanczosSolverType lanczosSolver(matrix,params);

			TridiagonalMatrixType ab;

			lanczosSolver.decomposition(modifVector,ab);
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
		PsimagLite::ProgressIndicator progress_;
		ParametersEngine<RealType> params_;
		RealType gsEnergy_;
		VectorType gsVector_; 
	}; // class ContinuedFraction
} // namespace Dmrg

#endif  // ENGINE_H_
