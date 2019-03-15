/*
*/

#ifndef LANCZOS_KITAEV_H
#define LANCZOS_KITAEV_H

#include "CrsMatrix.h"
#include "BitManip.h"
#include "TypeToString.h"
#include "SparseRow.h"
#include "ParametersKitaev.h"
#include "ModelBase.h"
#include "BasisKitaev.h"

namespace LanczosPlusPlus {

template<typename ComplexOrRealType,typename GeometryType,typename InputType>
class Kitaev  : public ModelBase<ComplexOrRealType,GeometryType,InputType> {

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef ModelBase<ComplexOrRealType,GeometryType,InputType> BaseType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	enum {SPIN_UP = ProgramGlobals::SPIN_UP, SPIN_DOWN = ProgramGlobals::SPIN_DOWN};

public:

	typedef ParametersKitaev<RealType,InputType> ParametersModelType;
	typedef BasisKitaev<GeometryType> BasisType;
	typedef typename BasisType::LabeledOperatorType LabeledOperatorType;
	typedef typename BasisType::BaseType BasisBaseType;
	typedef typename BasisType::WordType WordType;
	typedef typename BaseType::SparseMatrixType SparseMatrixType;
	typedef typename BaseType::VectorType VectorType;
	typedef PsimagLite::SparseRow<SparseMatrixType> SparseRowType;

	static const SizeType TWICE_THE_SPIN = BasisType::TWICE_THE_SPIN;

	Kitaev(InputType& io,
	       const GeometryType& geometry)
	    : mp_(io),
	      geometry_(geometry),
	      basis_(geometry),
	      jpp_(geometry_.numberOfSites(),geometry_.numberOfSites()),
	      jpm_(geometry_.numberOfSites(),geometry_.numberOfSites()),
	      jzz_(geometry_.numberOfSites(),geometry_.numberOfSites())
	{
		const RealType zeroPointTwentyFive = 0.25;
		SizeType n = geometry_.numberOfSites();

		if (geometry_.terms() != 3) {
			PsimagLite::String msg("Kitaev: must have 3 terms\n");
			throw PsimagLite::RuntimeError(msg);
		}

		for (SizeType i = 0; i < n; ++i) {
			for (SizeType j = 0; j < n; ++j) {
				// FIXME: MAKE SURE THAT WHEN i==j ALL CONNECTIONS ARE ZERO HERE
				ComplexOrRealType jxx = geometry_(i,0,j,0,0);
				ComplexOrRealType jyy = geometry_(i,0,j,0,1);
				jpm_(i, j) = zeroPointTwentyFive*(jxx + jyy);
				jpp_(i, j) = zeroPointTwentyFive*(jxx - jyy); // = jmm_(i,j)
				jzz_(i,j) = geometry_(i,0,j,0,2);
			}
		}
	}

	~Kitaev()
	{
		BaseType::deleteGarbage(garbage_);
	}

	SizeType size() const { return basis_.size(); }

	SizeType orbitals(SizeType) const { return 1; }

	void setupHamiltonian(SparseMatrixType& matrix) const
	{
		setupHamiltonian(matrix, basis_);
	}

	//! Gf. related functions below:
	void setupHamiltonian(SparseMatrixType& matrix,
	                      const BasisBaseType& basis) const
	{
		SizeType hilbert=basis.size();
		VectorRealType diag(hilbert,0.0);
		calcDiagonalElements(diag, basis);

		SizeType nsite = geometry_.numberOfSites();
		SizeType dummy = 0;
		SizeType orb = 0;

		matrix.resize(hilbert, hilbert);
		// Calculate off-diagonal elements AND store matrix
		SizeType nCounter=0;
		for (SizeType ispace = 0; ispace < hilbert; ++ispace) {
			SparseRowType sparseRow;
			matrix.setRow(ispace, nCounter);

			WordType ket = basis(ispace, dummy);
			// Save diagonal
			sparseRow.add(ispace, diag[ispace]);

			// s+s-
			for (SizeType i = 0; i < nsite; ++i) {
				SizeType val1 = basis.getN(ket, dummy, i, dummy, orb);
				if (val1 != TWICE_THE_SPIN) {
					setSplusSminus(sparseRow, ket, i, val1 + 1, basis);
					setSplusSplus(sparseRow, ket, i, val1 + 1, basis);
				}

				if (val1 != 0)
					setSminusSminus(sparseRow, ket, i, val1 - 1, basis);
			}

			nCounter += sparseRow.finalize(matrix);
		}

		matrix.setRow(hilbert,nCounter);
		matrix.checkValidity();
		assert(isHermitian(matrix));
	}

	bool hasNewParts(std::pair<SizeType,SizeType>& newParts,
	                 const std::pair<SizeType,SizeType>& oldParts,
	                 const LabeledOperatorType& lOperator,
	                 SizeType spin,
	                 SizeType orb) const
	{
		throw PsimagLite::RuntimeError("Kitaev::hasNewParts() unimplemented yet\n");
		//		if (what==ProgramGlobals::OPERATOR_SZ)
		//			return false;

		//		if (what == ProgramGlobals::OPERATOR_SPLUS ||
		//		        what == ProgramGlobals::OPERATOR_SMINUS)
		//			return hasNewPartsSplusOrMinus(newParts,oldParts,what,spin);

		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) +  "\n";
		str += PsimagLite::String("hasNewParts: unsupported operator ");
		str += lOperator.toString() + "\n";
		throw PsimagLite::RuntimeError(str.c_str());
	}

	const GeometryType& geometry() const { return geometry_; }

	const BasisType& basis() const { return basis_; }

	void printBasis(std::ostream& os) const { os<<basis_; }

	PsimagLite::String name() const { return __FILE__; }

	BasisType* createBasis(SizeType, SizeType) const
	{
		BasisType* ptr = new BasisType(geometry_);
		garbage_.push_back(ptr);
		return ptr;
	}

	void print(std::ostream& os) const { os<<mp_; }

	void printOperators(std::ostream& os) const
	{
		SizeType sites = geometry_.numberOfSites();
		for (SizeType site = 0; site < sites; ++site)
			printOperatorSz(site,os);
	}

private:

	void printOperatorSz(SizeType site, std::ostream& os) const
	{
		MatrixType matrix;
		setupOperator(matrix,"sz",site);
		os<<"#Operator_Sz_"<<"_"<<site<<"\n";
		os<<"#Matrix\n";
		os<<matrix;
	}

	void setupOperator(MatrixType& matrix,
	                   PsimagLite::String operatorName,
	                   SizeType site) const
	{
		SizeType hilbert = basis_.size();
		SizeType nsite = geometry_.numberOfSites();

		if (operatorName != "sz") {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "operator " + operatorName + " is unimplemented for this model\n";
			throw PsimagLite::RuntimeError(str);
		}

		if (site >= nsite) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "site requested " + ttos(site);
			str += " but number of sites= " + ttos(nsite) + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		matrix.resize(hilbert,hilbert);
		matrix.setTo(0.0);

		SizeType orb = 0;
		SizeType dummy = 0;
		for (SizeType ispace=0;ispace<hilbert;ispace++) {
			WordType ket = basis_(ispace,dummy);
			// assumes OPERATOR_SZ
			SizeType val1 = basis_.getN(ket,dummy,site,dummy,orb);
			RealType tmp = val1 - TWICE_THE_SPIN*0.5;

			matrix(ispace,ispace) = tmp;
		}
	}

	bool hasNewPartsSplusOrMinus(std::pair<SizeType,SizeType>& newParts,
	                             const std::pair<SizeType,SizeType>& oldParts,
	                             const LabeledOperatorType& lOperator,
	                             SizeType) const
	{
		if (TWICE_THE_SPIN != 1)
			throw PsimagLite::RuntimeError("BasisKitaev::hasNewPartsSplusOrMinus\n");

		newParts.first = oldParts.first;
		bool flag = false;
		if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SPLUS) {
			newParts.second = oldParts.second + 1;
			if (newParts.second > geometry_.numberOfSites())
				flag = true;
		} else if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SMINUS) {
			if (basis_.szPlusConst() == 0) flag =true;
			newParts.second = oldParts.second - 1;
		}

		if (!flag) return true;

		throw PsimagLite::RuntimeError("BasisKitaev:: S+ to all ups or S- to all downs\n");
	}

	void calcDiagonalElements(VectorRealType& diag,
	                          const BasisBaseType& basis) const
	{
		SizeType hilbert=basis.size();
		SizeType nsite = geometry_.numberOfSites();
		SizeType orb = 0;
		SizeType dummy = 0;

		// Calculate diagonal elements
		for (SizeType ispace = 0; ispace < hilbert; ++ispace) {
			WordType ket = basis(ispace, dummy);
			ComplexOrRealType s = 0;
			for (SizeType i = 0; i < nsite; ++i) {
				SizeType val1 = basis.getN(ket,dummy,i,dummy,orb);
				RealType tmp1 = val1 - TWICE_THE_SPIN*0.5;
				if (i < mp_.magneticField.size()) s += mp_.magneticField[i]*tmp1;

				for (SizeType j = i + 1; j < nsite; ++j) {
					SizeType val2 = basis.getN(ket, dummy, j, dummy, orb);
					RealType tmp2 = val2 - TWICE_THE_SPIN*0.5;

					// Sz Sz term:
					s += tmp1*tmp2*jzz_(i,j);
				}
			}

			assert(fabs(PsimagLite::imag(s))<1e-12);
			diag[ispace] = PsimagLite::real(s);
		}
	}

	void setSplusSminus(SparseRowType &sparseRow,
	                    const WordType& ket,
	                    SizeType i,
	                    SizeType val1,
	                    const BasisBaseType& basis) const
	{
		const RealType zero = 0;
		SizeType nsite = geometry_.numberOfSites();
		SizeType dummy = 0;
		SizeType orb = 0;
		RealType spin = TWICE_THE_SPIN*0.5;

		for (SizeType j = 0; j < nsite; ++j) {
			if (i == j) continue;
			if (jpm_(i,j) == zero) continue;

			SizeType val2 = basis.getN(ket, dummy, j, dummy, orb);
			if (val2 == 0) continue;
			RealType m2 = val2 - spin;
			--val2;
			RealType m1 = val2 - spin;
			WordType bra;

			basis.getBra(bra,ket,i,LabeledOperatorType(val1),j,val2);
			SizeType temp = basis.perfectIndex(bra,dummy);
			RealType tmp = sqrt(spin*(spin+1.0) - m1*(m1+1.0));
			tmp *= sqrt(spin*(spin+1.0) - m2*(m2-1.0));
			sparseRow.add(temp,tmp*jpm_(i,j));
		}
	}

	void setSplusSplus(SparseRowType &sparseRow,
	                   const WordType& ket,
	                   SizeType i,
	                   SizeType val1,
	                   const BasisBaseType& basis) const
	{
		const RealType zero = 0;
		SizeType nsite = geometry_.numberOfSites();
		SizeType dummy = 0;
		SizeType orb = 0;
		if (TWICE_THE_SPIN != 1)
			err("Kitaev only for spin=1/2 for now (sorry)\n");

		for (SizeType j = 0; j < nsite; ++j) {
			if (i == j) continue;
			if (jpp_(i, j) == zero) continue;

			SizeType val2 = basis.getN(ket, dummy, j, dummy, orb);

			if (val2 == TWICE_THE_SPIN) continue;

			WordType bra;
			basis.getBra(bra, ket, i, val1, j, val2 + 1);
			SizeType temp = basis.perfectIndex(bra, dummy);
			RealType tmp = 1.0; // FIXME: CORRECT FOR Kitaev with spin != 1/2
			sparseRow.add(temp, tmp*jpp_(i, j));
		}
	}

	void setSminusSminus(SparseRowType &sparseRow,
	                    const WordType& ket,
	                    SizeType i,
	                    SizeType val1,
	                    const BasisBaseType& basis) const
	{
		if (TWICE_THE_SPIN != 1)
			err("Kitaev only for spin=1/2 for now (sorry)\n");
		setSplusSminus(sparseRow, ket, i, val1, basis);
	}

	RealType signSplusSminus(SizeType i,
	                         SizeType j,
	                         const WordType& bra1,
	                         const WordType& bra2) const
	{
		SizeType n = geometry_.numberOfSites();
		int s = 1;
		if (j>0) s *= parityFrom(0,j-1,bra2);
		if (i>0) s *= parityFrom(0,i-1,bra2);
		if (i<n-1) s *= parityFrom(i+1,n-1,bra1);
		if (j<n-1) s *= parityFrom(j+1,n-1,bra1);
		return s;
	}

	// from i to j including i and j
	// assumes i<=j
	int parityFrom(SizeType i,SizeType j,const WordType& ket) const
	{
		if (i==j) return (BasisType::bitmask(j) & ket) ? -1 : 1;
		assert(i<j);
		//j>i>=0 now
		WordType mask = ket;
		mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
		int s=(PsimagLite::BitManip::count(mask) & 1) ? -1 : 1;
		// Is there something of this species at i?
		if (BasisType::bitmask(i) & ket) s = -s;
		// Is there something of this species at j?
		if (BasisType::bitmask(j) & ket) s = -s;
		return s;
	}

	const ParametersModelType mp_;
	const GeometryType& geometry_;
	BasisType basis_;
	PsimagLite::Matrix<ComplexOrRealType> jpp_;
	PsimagLite::Matrix<ComplexOrRealType> jpm_;
	PsimagLite::Matrix<ComplexOrRealType> jzz_;
	mutable typename PsimagLite::Vector<BasisType*>::Type garbage_;
}; // class Kitaev
} // namespace LanczosPlusPlus
#endif

