/*
*/

#ifndef HUBBARDLANCZOS_H
#define HUBBARDLANCZOS_H

#include "BasisHubbardLanczos.h"
#include "BitManip.h"
#include "TypeToString.h"
#include "SparseRow.h"
#include "ParametersModelHubbard.h"
#include "ProgramGlobals.h"
#include "../../Engine/ModelBase.h"

namespace LanczosPlusPlus {

template<typename ComplexOrRealType,typename GeometryType,typename InputType>
class HubbardOneOrbital : public ModelBase<ComplexOrRealType,GeometryType,InputType> {

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef ModelBase<ComplexOrRealType,GeometryType,InputType> BaseType;

	enum {TERM_HOPPING=0,TERM_NINJ=1,TERM_SUPER=2};

	enum {SO_HOPPING_TERM = 1};

	enum {SPIN_UP = ProgramGlobals::SPIN_UP, SPIN_DOWN = ProgramGlobals::SPIN_DOWN};

public:

	typedef ParametersModelHubbard<RealType,InputType> ParametersModelType;
	typedef BasisHubbardLanczos<GeometryType> BasisType;
	typedef typename BasisType::PairIntType PairIntType;
	typedef typename BasisType::BaseType BasisBaseType;
	typedef typename BasisType::WordType WordType;
	typedef typename BaseType::VectorSizeType VectorSizeType;
	typedef typename BaseType::SparseMatrixType SparseMatrixType;
	typedef typename BaseType::VectorType VectorType;
	typedef PsimagLite::SparseRow<SparseMatrixType> SparseRowType;

	static int const FERMION_SIGN = BasisType::FERMION_SIGN;

	HubbardOneOrbital(SizeType nup,
	                  SizeType ndown,
	                  const ParametersModelType& mp,
	                  const GeometryType& geometry)
	    : mp_(mp),
	      geometry_(geometry),
	      basis_(geometry,nup,ndown),
	      hoppings_(geometry_.numberOfSites(),geometry_.numberOfSites()),
	      hasJcoupling_(false),
	      hasCoulombCoupling_(false)
	{
		if (mp_.model == "HubbardOneBandExtended" ||
		    mp_.model == "SuperHubbardExtended") hasCoulombCoupling_ = true;

		if (mp_.model == "SuperHubbardExtended") hasJcoupling_ = true;

		bool hasSpinOrbitKaneMele = (mp_.model == "KaneMeleHubbard");

		if (hasCoulombCoupling_ && geometry_.terms()<2) {
			throw PsimagLite::RuntimeError("HubbardOneOrbital::ctor(): ColoumbCoupling\n");
		}

		if (hasJcoupling_ && geometry_.terms()<3) {
			throw PsimagLite::RuntimeError("HubbardOneOrbital::ctor(): jCoupling\n");
		}

		if (hasSpinOrbitKaneMele && geometry_.terms() != 2) {
			throw PsimagLite::RuntimeError("HubbardOneOrbital::ctor(): KaneMeleHubbard");
		}

		SizeType n = geometry_.numberOfSites();
		for (SizeType j=0;j<n;j++) {
			for (SizeType i=0;i<j;i++) {

				hoppings_(i,j) = geometry_(i,0,j,0,TERM_HOPPING);

				if (hasSpinOrbitKaneMele)
					hoppings_(i,j) += geometry_(i,0,j,0,SO_HOPPING_TERM);
			}
		}

		for (SizeType j=0;j<n;j++) {
			for (SizeType i=j+1;i<n;i++) {
				hoppings_(i,j) = std::conj(hoppings_(i,j));
			}
		}
	}

	~HubbardOneOrbital()
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
		setupHamiltonian(matrix,basis_);
	}

	//! Gf. related functions below:
	void setupHamiltonian(SparseMatrixType& matrix,
	                      const BasisBaseType& basis) const
	{
		SizeType hilbert=basis.size();
		typename PsimagLite::Vector<RealType>::Type diag(hilbert);
		calcDiagonalElements(diag,basis);

		SizeType nsite = geometry_.numberOfSites();

		matrix.resize(hilbert,hilbert);
		// Calculate off-diagonal elements AND store matrix
		SizeType nCounter=0;
		for (SizeType ispace=0;ispace<hilbert;ispace++) {
			SparseRowType sparseRow;
			matrix.setRow(ispace,nCounter);
			WordType ket1 = basis(ispace,SPIN_UP);
			WordType ket2 = basis(ispace,SPIN_DOWN);
			// Save diagonal
			sparseRow.add(ispace,diag[ispace]);
			for (SizeType i=0;i<nsite;i++) {
				setHoppingTerm(sparseRow,ket1,ket2,i,basis);
				setJTermOffDiagonal(sparseRow,ket1,ket2,i,basis);
			}

			nCounter += sparseRow.finalize(matrix);
		}

		matrix.setRow(hilbert,nCounter);
	}

	bool hasNewParts(std::pair<SizeType,SizeType>& newParts,
	                 SizeType what,
	                 SizeType spin,
	                 const std::pair<SizeType,SizeType>& orbs) const
	{
		if (what==ProgramGlobals::OPERATOR_C || what==ProgramGlobals::OPERATOR_CDAGGER)
			return hasNewPartsCorCdagger(newParts,what,spin,orbs);
		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) +  "\n";
		str += PsimagLite::String("hasNewParts: unsupported operator ");
		str += ProgramGlobals::id2Operator(what) + "\n";
		throw std::runtime_error(str.c_str());
	}

	const GeometryType& geometry() const { return geometry_; }

	const BasisType& basis() const { return basis_; }

	PsimagLite::String name() const { return __FILE__; }

	BasisType* createBasis(SizeType nup, SizeType ndown) const
	{
		BasisType* ptr = new BasisType(geometry_,nup,ndown);
		garbage_.push_back(ptr);
		return ptr;
	}

	void print(std::ostream& os) const { os<<mp_; }

	void printOperators(std::ostream& os) const
	{
		SizeType spin = SPIN_UP;
		SizeType site = 0;
		SizeType nup = basis_.electrons(SPIN_UP);
		SizeType ndown = basis_.electrons(SPIN_DOWN);
		if (nup == 0) {
			os<<"#Operator_c_"<<spin<<"_"<<site<<"\n";
			os<<"0 0\n";
			return;
		}

		BasisType*  basis = createBasis(nup-1, ndown);
		VectorSizeType opt(2,0);
		opt[0] = site;
		opt[1] = spin;
		SparseMatrixType matrix;
		setupOperator(matrix,*basis,"c",opt);
		MatrixType fm;
		crsMatrixToFullMatrix(fm,matrix);
		os<<"#Operator_c_"<<spin<<"_"<<site<<"\n";
		os<<fm;
	}

private:

	void setupOperator(SparseMatrixType& matrix,
	                   const BasisBaseType& basis,
	                   PsimagLite::String operatorName,
	                   const VectorSizeType& operatorOptions) const
	{
		SizeType hilbert=basis.size();
		SizeType nsite = geometry_.numberOfSites();
		if (operatorOptions.size() < 2) {
			if (operatorName != "i") {
				PsimagLite::String str(__FILE__);
				str += " " + ttos(__LINE__) + "\n";
				str += "operator "+ operatorName + " needs at least two options\n";
				throw PsimagLite::RuntimeError(str);
			}

			matrix.makeDiagonal(hilbert,1.0);
			return;
		}

		SizeType site = operatorOptions[0];
		if (site >= nsite) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "site requested " + ttos(site);
			str += " but number of sites= " + ttos(nsite) + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		SizeType id = 0;
		if (operatorName == "c") {
			id = ProgramGlobals::OPERATOR_C;
		} else {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "operator " + operatorName + " is unimplemented for this model\n";
			throw PsimagLite::RuntimeError(str);
		}

		SizeType spin = operatorOptions[1];
		matrix.resize(hilbert,hilbert);
		// Calculate off-diagonal elements AND store matrix
		SizeType nCounter=0;
		for (SizeType ispace=0;ispace<hilbert;ispace++) {
			SparseRowType sparseRow;
			matrix.setRow(ispace,nCounter);
			WordType ket1 = basis_(ispace,SPIN_UP);
			WordType ket2 = basis_(ispace,SPIN_DOWN);
			// assumes OPERATOR_C
			PairIntType bra = basis.getBraIndex(ket1,ket2,id,site,spin,0);
			if (bra.first < 0) continue;
			sparseRow.add(bra.first,bra.second);
			nCounter += sparseRow.finalize(matrix);
		}

		matrix.setRow(hilbert,nCounter);
		matrix.checkValidity();
	}

	bool hasNewPartsCorCdagger(std::pair<SizeType,SizeType>& newParts,
	                           SizeType what,
	                           SizeType spin,
	                           const std::pair<SizeType,SizeType>&) const
	{
		int newPart1=basis_.electrons(SPIN_UP);
		int newPart2=basis_.electrons(SPIN_DOWN);
		int c = (what==ProgramGlobals::OPERATOR_C) ? -1 : 1;
		if (spin==SPIN_UP) newPart1 += c;
		else newPart2 += c;

		if (newPart1<0 || newPart2<0) return false;
		SizeType nsite = geometry_.numberOfSites();
		if (SizeType(newPart1)>nsite || SizeType(newPart2)>nsite) return false;
		if (newPart1==0 && newPart2==0) return false;
		newParts.first = SizeType(newPart1);
		newParts.second = SizeType(newPart2);
		return true;
	}

	void calcDiagonalElements(typename PsimagLite::Vector<RealType>::Type& diag,
	                          const BasisBaseType& basis) const
	{
		SizeType hilbert=basis.size();
		SizeType nsite = geometry_.numberOfSites();
		SizeType orb = 0;

		// Calculate diagonal elements
		for (SizeType ispace=0;ispace<hilbert;ispace++) {
			WordType ket1 = basis(ispace,SPIN_UP);
			WordType ket2 = basis(ispace,SPIN_DOWN);
			ComplexOrRealType s=0;
			for (SizeType i=0;i<nsite;i++) {

				// Hubbard term U0
				s += mp_.hubbardU[i] *
				        basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb) *
				        basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb);

				// SzSz
				for (SizeType j=0;j<nsite;j++) {
					ComplexOrRealType value = jCoupling(i,j);
					if (std::real(value) == 0 && std::imag(value) == 0) continue;
					s += value*0.5* // double counting i,j
					        szTerm(ket1,ket2,i,basis)*
					        szTerm(ket1,ket2,j,basis);
				}

				// Coulomb
				RealType ne = (basis.getN(ket1,ket2,i,SPIN_UP,orb) +
				               basis.getN(ket1,ket2,i,SPIN_DOWN,orb));

				for (SizeType j=0;j<nsite;j++) {
					ComplexOrRealType value = coulombCoupling(i,j);
					if (std::real(value) == 0 && std::imag(value) == 0) continue;
					RealType tmp2 = basis.getN(ket1,ket2,j,SPIN_UP,orb) +
					        basis.getN(ket1,ket2,j,SPIN_DOWN,orb);
					s += value * ne * tmp2;
				}

				// Potential term
				RealType tmp = mp_.potentialV[i];
				if (mp_.potentialT.size()>0)
					tmp += mp_.potentialT[i]*mp_.timeFactor;
				if (tmp!=0) s += tmp * ne;
			}

			assert(fabs(std::imag(s))<1e-12);
			diag[ispace] = std::real(s);
		}
	}

	void setHoppingTerm(SparseRowType& sparseRow,
	                    const WordType& ket1,
	                    const WordType& ket2,
	                    SizeType i,
	                    const BasisBaseType& basis) const
	{
		WordType s1i=(ket1 & BasisType::bitmask(i));
		if (s1i>0) s1i=1;
		WordType s2i=(ket2 & BasisType::bitmask(i));
		if (s2i>0) s2i=1;

		SizeType nsite = geometry_.numberOfSites();
		SizeType orb = 0;

		// Hopping term
		for (SizeType j=0;j<nsite;j++) {
			if (j<i) continue;
			ComplexOrRealType h = hoppings_(i,j);
			if (std::real(h) == 0 && std::imag(h) == 0) continue;
			WordType s1j= (ket1 & BasisType::bitmask(j));
			if (s1j>0) s1j=1;
			WordType s2j= (ket2 & BasisType::bitmask(j));
			if (s2j>0) s2j=1;

			if (s1i+s1j==1) {
				WordType bra1= ket1 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
				SizeType temp = basis.perfectIndex(bra1,ket2);
				RealType extraSign = (s1i==1) ? FERMION_SIGN : 1;
				RealType tmp2 = basis.doSign(ket1,ket2,i,orb,j,orb,SPIN_UP);
				ComplexOrRealType cTemp = h*extraSign*tmp2;
				if (s1i == 0) cTemp = std::conj(cTemp);
				assert(temp<basis.size());
				sparseRow.add(temp,cTemp);
			}

			if (s2i+s2j==1) {
				WordType bra2= ket2 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
				SizeType temp = basis.perfectIndex(ket1,bra2);
				RealType extraSign = (s2i==1) ? FERMION_SIGN : 1;
				RealType tmp2 = basis.doSign(ket1,ket2,i,orb,j,orb,SPIN_DOWN);
				ComplexOrRealType cTemp = h*extraSign*tmp2;
				if (s2j == 0) cTemp = std::conj(cTemp);
				assert(temp<basis.size());
				sparseRow.add(temp,cTemp);
			}
		}
	}

	void setJTermOffDiagonal(SparseRowType& sparseRow,
	                         const WordType& ket1,
	                         const WordType& ket2,
	                         SizeType i,
	                         const BasisBaseType& basis) const
	{
		for (SizeType j=0;j<geometry_.numberOfSites();j++) {
			ComplexOrRealType value = jCoupling(i,j)*0.5;
			if (std::real(value) == 0 && std::imag(value) == 0) continue;
			value *= 0.5; // double counting i,j
			assert(i!=j);

			RealType sign = jTermSign(ket1,ket2,i,j,basis);
			setSplusSminus(sparseRow,ket1,ket2,i,j,value*sign,basis);
			setSplusSminus(sparseRow,ket1,ket2,j,i,value*sign,basis);
		}
	}

	void setSplusSminus(SparseRowType &sparseRow,
	                    const WordType& ket1,
	                    const WordType& ket2,
	                    SizeType i,
	                    SizeType j,
	                    ComplexOrRealType value,
	                    const BasisBaseType &basis) const
	{
		if (splusSminusNonZero(ket1,ket2,i,j,basis)==0) return;

		assert(i!=j);
		WordType bra1 = ket1 ^ (BasisType::bitmask(i)|BasisType::bitmask(j));
		WordType bra2 = ket2 ^ (BasisType::bitmask(i)|BasisType::bitmask(j));
		SizeType temp = basis.perfectIndex(bra1,bra2);
		sparseRow.add(temp,value);
	}

	SizeType splusSminusNonZero(const WordType& ket1,
	                            const WordType& ket2,
	                            SizeType i,
	                            SizeType j,
	                            const BasisBaseType &basis) const
	{
		SizeType orb = 0;
		if (basis.isThereAnElectronAt(ket1,ket2,j,SPIN_UP,orb)==0) return 0;
		if (basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb)==1) return 0;
		if (basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb)==0) return 0;
		if (basis.isThereAnElectronAt(ket1,ket2,j,SPIN_DOWN,orb)==1) return 0;
		return 1;
	}

	int jTermSign(const WordType& ket1,
	              const WordType& ket2,
	              SizeType i,
	              SizeType j,
	              const BasisBaseType& basis) const
	{
		if (i>j) return jTermSign(ket1,ket2,j,i,basis);
		SizeType orb = 0;
		int x = basis.doSign(ket1,ket2,i,orb,j,orb,SPIN_UP);
		x *= basis.doSign(ket1,ket2,i,orb,j,orb,SPIN_DOWN);
		return x;
	}

	RealType szTerm(const WordType& ket1,
	                const WordType& ket2,
	                SizeType i,
	                const BasisBaseType& basis) const
	{
		SizeType orb = 0;
		RealType sz = basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb);
		sz -= basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb);
		return 0.5*sz;
	}

	ComplexOrRealType jCoupling(SizeType i,SizeType j) const
	{
		if (!hasJcoupling_) return 0;
		return geometry_(i,0,j,0,TERM_SUPER);
	}

	ComplexOrRealType coulombCoupling(SizeType i,SizeType j) const
	{
		if (!hasCoulombCoupling_) return 0;
		return geometry_(i,0,j,0,TERM_NINJ);
	}

	const ParametersModelType mp_;
	const GeometryType& geometry_;
	BasisType basis_;
	PsimagLite::Matrix<ComplexOrRealType> hoppings_;
	bool hasJcoupling_;
	bool hasCoulombCoupling_;
	mutable typename PsimagLite::Vector<BasisType*>::Type garbage_;
}; // class HubbardOneOrbital
} // namespace LanczosPlusPlus
#endif

