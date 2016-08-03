
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
#include "DefaultSymmetry.h"
#include "TypeToString.h"

namespace LanczosPlusPlus {
template<typename ModelType_,
         template<typename,typename> class InternalProductTemplate,
         typename SpecialSymmetryType>
class Engine  {

public:

	enum {SPIN_UP = ProgramGlobals::SPIN_UP, SPIN_DOWN = ProgramGlobals::SPIN_DOWN};

	enum {OPERATOR_N = ProgramGlobals::OPERATOR_N};

	typedef ModelType_ ModelType;
	typedef typename ModelType::InputType InputType;
	typedef typename ModelType::SparseMatrixType SparseMatrixType;
	typedef typename ModelType::BasisBaseType BasisType;
	typedef InternalProductTemplate<ModelType,SpecialSymmetryType> InternalProductType;
	typedef DefaultSymmetry<typename ModelType::GeometryType,BasisType> DefaultSymmetryType;
	typedef InternalProductTemplate<ModelType,DefaultSymmetryType> InternalProductDefaultType;
	typedef typename SpecialSymmetryType::GeometryType GeometryType;
	typedef typename GeometryType::ComplexOrRealType ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Random48<RealType> RandomType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef PsimagLite::LanczosSolver<ParametersForSolverType,
	                                  InternalProductType,
	                                  VectorType> LanczosSolverType;
	typedef PsimagLite::LanczosSolver<ParametersForSolverType,
	                                  InternalProductDefaultType,
	                                  VectorType> LanczosSolverDefaultType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	// ContF needs to support concurrency FIXME
	static const SizeType parallelRank_ = 0;
	static const SizeType CHECK_HERMICITY = 1;

	enum {PLUS,MINUS};

	Engine(const ModelType& model,
	       SizeType,
	       InputType& io)
	    : model_(model),
	      progress_("Engine"),
	      io_(io),
	      options_("")
	{
		io_.readline(options_,"SolverOptions=");
		computeGroundState();
	}

	RealType gsEnergy() const
	{
		return gsEnergy_;
	}

	const VectorType& eigenvector() const
	{
		return gsVector_;
	}

	//! Calc Green function G(isite,jsite)  (still diagonal in spin)
	template<typename ContinuedFractionCollectionType>
	void spectralFunction(ContinuedFractionCollectionType& cfCollection,
	                      VectorStringType& vstr,
	                      SizeType what2,
	                      int isite,
	                      int jsite,
	                      const PsimagLite::Vector<PairType>::Type& spins,
	                      const PairType& orbs) const
	{
		std::cout<<"orbitals="<<orbs.first<<" "<<orbs.second<<"\n";
		for (SizeType i=0;i<spins.size();i++) {
			std::cout<<"spins="<<spins[i].first<<" "<<spins[i].second<<"\n";
			spectralFunction(cfCollection,vstr,what2,isite,jsite,spins[i],orbs);
		}
	}

	/* PSIDOC SpectralFunctions
	Here we document the spectral functions and Green function G(isite,jsite)  
	(still diagonal in spin)
	*/
	template<typename ContinuedFractionCollectionType>
	void spectralFunction(ContinuedFractionCollectionType& cfCollection,
	                      VectorStringType& vstr,
	                      SizeType what2,
	                      int isite,
	                      int jsite,
	                      const PairType& spins,
	                      const PairType& orbs) const
	{
		if (spins.first!=spins.second) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "spectralFunction: no support yet for off-diagonal spin\n";
			throw std::runtime_error(str.c_str());
		}

		typedef typename ContinuedFractionCollectionType::ContinuedFractionType
		        ContinuedFractionType;

		const BasisType* basisNew = 0;
		bool isDiagonal = (isite==jsite && orbs.first==orbs.second);

		for (SizeType type=0;type<4;type++) {
			if (isDiagonal && type>1) continue;

			SizeType operatorLabel= (type&1) ?  what2 : ProgramGlobals::transposeConjugate(what2);
			if (ProgramGlobals::needsNewBasis(operatorLabel)) {
				assert(spins.first==spins.second);
				std::pair<SizeType,SizeType> newParts(0,0);
				if (!model_.hasNewParts(newParts,operatorLabel,spins.first,orbs)) continue;
				// Create new bases
				basisNew = model_.createBasis(newParts.first,newParts.second);
			} else {
				basisNew = &model_.basis();
			}
			VectorType modifVector;
			getModifiedState(modifVector,
			                 operatorLabel,
			                 gsVector_,
			                 *basisNew,
			                 type,
			                 isite,
			                 jsite,
			                 spins.first,
			                 orbs);

			DefaultSymmetryType symm(*basisNew,model_.geometry(),"");
			InternalProductTemplate<ModelType,DefaultSymmetryType> matrix(model_,
			                                                              *basisNew,
			                                                              symm);
			ContinuedFractionType cf(cfCollection.freqType());

			if (PsimagLite::norm(modifVector)<1e-10) {
				std::cerr<<"spectralFunction: modifVector==0, type="<<type<<"\n";
			}

			calcSpectral(cf,operatorLabel,modifVector,matrix,type,spins.first,isDiagonal);
			PsimagLite::String str = ttos(spins.first) + "," + ttos(type) + ",";
			str += ttos(orbs.first) + "," + ttos(orbs.second);
			vstr.push_back(str);
			cfCollection.push(cf);
		}
	}

	void twoPoint(PsimagLite::Matrix<typename VectorType::value_type>& result,
	              SizeType what2,
	              const PsimagLite::Vector<PairType>::Type& spins,
	              const PairType& orbs) const
	{
		for (SizeType i=0;i<spins.size();i++) {
			std::cout<<"spins="<<spins[i].first<<" "<<spins[i].second<<"\n";
			twoPoint(result,what2,spins[i],orbs);
		}

	}

	/* PSIDOC TwoPointCorrelations
	Here we document the two-point correlations
	*/
	void twoPoint(PsimagLite::Matrix<typename VectorType::value_type>& result,
	              SizeType what2,
	              const PairType& spins,
	              const PairType& orbs) const
	{
		const BasisType* basisNew = 0;

		if (ProgramGlobals::needsNewBasis(what2)) {
			if (spins.first!=spins.second) {
				PsimagLite::String str(__FILE__);
				str += " " + ttos(__LINE__) + "\n";
				str += "twoPoint: no support yet for off-diagonal spin ";
				str += "when needs new basis\n";
				throw std::runtime_error(str.c_str());
			}

			std::pair<SizeType,SizeType> newParts(0,0);
			if (!model_.hasNewParts(newParts,what2,spins.first,orbs)) return;

			basisNew = model_.createBasis(newParts.first,newParts.second);

			std::cerr<<"basisNew.size="<<basisNew->size()<<" ";
			std::cerr<<"newparts.first="<<newParts.first<<" ";
			std::cerr<<"newparts.second="<<newParts.second<<"\n";
		} else {
			basisNew = &model_.basis();
		}

		SizeType total =result.n_row();

		for (SizeType isite=0;isite<total;isite++)
			for (SizeType jsite=0;jsite<total;jsite++)
				result(isite,jsite) = -100;

		RealType isign = 1.0;

		typename VectorType::value_type sum = 0;
		std::cout<<"orbs="<<orbs.first<<" "<<orbs.second<<"\n";
		for (SizeType isite=0;isite<total;isite++) {
			VectorType modifVector1(basisNew->size(),0);
			if (orbs.first>=model_.orbitals(isite)) continue;
			accModifiedState(modifVector1,what2,*basisNew,gsVector_,
			                 isite,spins.first,orbs.first,isign);
			for (SizeType jsite=0;jsite<total;jsite++) {
				VectorType modifVector2(basisNew->size(),0);
				if (orbs.second>=model_.orbitals(jsite)) continue;
				accModifiedState(modifVector2,what2,*basisNew,gsVector_,
				                 jsite,spins.second,orbs.second,isign);
				result(isite,jsite) =  modifVector2*modifVector1;
				if (isite==jsite) sum += result(isite,isite);
			}
		}
		std::cout<<"MatrixDiagonal = "<<sum<<"\n";
	}

private:

	void accModifiedState_(VectorType &z,
	                       SizeType operatorLabel,
	                       const BasisType& newBasis,
	                       const VectorType& gsVector,
	                       SizeType site,
	                       SizeType spin,
	                       SizeType orb,
	                       RealType isign) const
	{
		for (SizeType ispace=0;ispace<model_.basis().size();ispace++) {
			ProgramGlobals::WordType ket1 = model_.basis()(ispace,SPIN_UP);
			ProgramGlobals::WordType ket2 = model_.basis()(ispace,SPIN_DOWN);
			ProgramGlobals::PairIntType tempValue = newBasis.getBraIndex(ket1,
			                                                             ket2,
			                                                             operatorLabel,
			                                                             site,
			                                                             spin,
			                                                             orb);
			int temp = tempValue.first;
			int value = tempValue.second;
			if (temp>=0 && SizeType(temp)>=z.size()) {
				PsimagLite::String s = "old basis=" + ttos(model_.basis().size());
				s += " newbasis=" + ttos(newBasis.size());
				s += "\n";
				s += "operatorLabel=" + ttos(operatorLabel) + " spin=" + ttos(spin);
				s += " site=" + ttos(site);
				s += "ket1=" + ttos(ket1) + " and ket2=" + ttos(ket2);
				s += "\n";
				s += "getModifiedState: z.size=" + ttos(z.size());
				s += " but temp=" + ttos(temp) + "\n";
				throw std::runtime_error(s.c_str());
			}
			if (temp<0) continue;
			int mysign = (ProgramGlobals::isFermionic(operatorLabel)) ?
			            model_.basis().doSignGf(ket1,ket2,site,spin,orb) : 1;
			if (operatorLabel == ProgramGlobals::OPERATOR_SPLUS ||
			        operatorLabel == ProgramGlobals::OPERATOR_SMINUS)
				mysign *= model_.basis().doSignSpSm(ket1,ket2,site,spin,orb);

			z[temp] += isign*mysign*value*gsVector[ispace];
		}
	}

	void getModifiedState(VectorType& modifVector,
	                      SizeType operatorLabel,
	                      const VectorType& gsVector,
	                      const BasisType& basisNew,
	                      SizeType type,
	                      SizeType isite,
	                      SizeType jsite,
	                      SizeType spin,
	                      const std::pair<SizeType,SizeType>& orbs) const
	{
		modifVector.resize(basisNew.size());
		for (SizeType temp=0;temp<modifVector.size();temp++)
			modifVector[temp]=0.0;

		accModifiedState_(modifVector,
		                  operatorLabel,
		                  basisNew,
		                  gsVector,
		                  isite,
		                  spin,
		                  orbs.first,
		                  1.0);
		std::cerr<<"isite="<<isite<<" type="<<type;
		std::cerr<<" modif="<<(modifVector*modifVector)<<"\n";
		if (model_.name()=="Tj1Orb.h" && isite==jsite) return;

		RealType isign= (type > 1) ? -1.0 : 1.0;
		accModifiedState_(modifVector,
		                  operatorLabel,
		                  basisNew,
		                  gsVector,
		                  jsite,
		                  spin,
		                  orbs.second,
		                  isign);
		std::cerr<<"jsite="<<jsite<<" type="<<type;
		std::cerr<<" modif="<<(modifVector*modifVector)<<"\n";
	}

	void accModifiedState(VectorType& z,
	                      SizeType operatorLabel,
	                      const BasisType& newBasis,
	                      const VectorType& gsVector,
	                      SizeType site,
	                      SizeType spin,
	                      SizeType orb,
	                      RealType isign) const
	{
		if (model_.name()=="Tj1Orb.h")
			accModifiedState_(z,operatorLabel,newBasis,gsVector,site,spin,orb,isign);

		if (operatorLabel==OPERATOR_N) {
			accModifiedState_(z,operatorLabel,newBasis,gsVector,site,spin,orb,isign);
			return;
		} else if (operatorLabel==ProgramGlobals::OPERATOR_SZ) {
			accModifiedState_(z,OPERATOR_N,newBasis,gsVector,site,SPIN_UP,orb,isign*0.5);
			accModifiedState_(z,OPERATOR_N,newBasis,gsVector,site,SPIN_DOWN,orb,-isign*0.5);
			return;
		}

		accModifiedState_(z,operatorLabel,newBasis,gsVector,site,spin,orb,isign);
	}

	void computeGroundState()
	{
		SpecialSymmetryType rs(model_.basis(),model_.geometry(),options_);
		InternalProductType hamiltonian(model_,rs);
		ParametersForSolverType params(io_,"Lanczos");
		LanczosSolverType lanczosSolver(hamiltonian,params);

		gsEnergy_ = 1e10;
		SizeType offset = model_.size();
		SizeType currentOffset = 0;
		for (SizeType i=0;i<rs.sectors();i++) {
			hamiltonian.specialSymmetrySector(i);
			VectorType gsVector1(hamiltonian.rank());
			if (gsVector1.size()==0) continue;
			RealType gsEnergy1 = 0;

			try {
				lanczosSolver.computeGroundState(gsEnergy1,gsVector1);
			} catch (std::exception& e) {

				std::cerr<<"Engine: Lanczos Solver failed ";
				std::cerr<<" trying exact diagonalization...\n";
				VectorRealType eigs(hamiltonian.rank());
				MatrixType fm;
				hamiltonian.fullDiag(eigs,fm);
				for (SizeType j = 0; j < eigs.size(); ++j)
					gsVector1[j] = fm(j,0);
				gsEnergy1 = eigs[0];
				std::cout<<"Found lowest eigenvalue= "<<gsEnergy1<<"\n";
			}

			if (gsEnergy1<gsEnergy_) {
				gsVector_=gsVector1;
				gsEnergy_=gsEnergy1;
				offset = currentOffset;
			}
			currentOffset +=  gsVector1.size();
		}
		rs.transformGs(gsVector_,offset);
		std::cout<<"#GSNorm="<<PsimagLite::real(gsVector_*gsVector_)<<"\n";
	}

	template<typename ContinuedFractionType>
	void calcSpectral(ContinuedFractionType& cf,
	                  SizeType what2,
	                  const VectorType& modifVector,
	                  const InternalProductDefaultType& matrix,
	                  SizeType type,
	                  SizeType,
	                  bool isDiagonal) const
	{
		typedef typename ContinuedFractionType::TridiagonalMatrixType
		        TridiagonalMatrixType;

		ParametersForSolverType params(io_,"Spectral");

		LanczosSolverDefaultType lanczosSolver(matrix,params);

		TridiagonalMatrixType ab;

		lanczosSolver.decomposition(modifVector,ab);
		typename VectorType::value_type weight = modifVector*modifVector;

		int s = (type&1) ? -1 : 1;
		RealType s2 = (type>1) ? -1 : 1;
		if (!ProgramGlobals::isFermionic(what2)) s2 *= s;
		RealType diagonalFactor = (isDiagonal) ? 1 : 0.5;
		s2 *= diagonalFactor;

		const MatrixRealType& reortho = lanczosSolver.reorthogonalizationMatrix();

		cf.set(ab,reortho,gsEnergy_,PsimagLite::real(weight*s2),s);

	}

	const ModelType& model_;
	PsimagLite::ProgressIndicator progress_;
	InputType& io_;
	PsimagLite::String options_;
	RealType gsEnergy_;
	VectorType gsVector_;
}; // class ContinuedFraction
} // namespace Dmrg

#endif  // ENGINE_H_

