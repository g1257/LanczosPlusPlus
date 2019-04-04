#ifndef HUBBARDHELPER_H
#define HUBBARDHELPER_H
#include "CrsMatrix.h"
#include "SparseRow.h"
#include "ProgramGlobals.h"

namespace LanczosPlusPlus {

template<typename ModelBaseType, typename BasisType, typename ModelParamsType>
class HubbardHelper {

public:

	typedef typename ModelBaseType::SparseMatrixType SparseMatrixType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::RealType RealType;
	typedef typename ModelBaseType::BasisBaseType BasisBaseType;
	typedef PsimagLite::SparseRow<SparseMatrixType> SparseRowType;
	typedef typename BasisType::WordType WordType;
	typedef typename ModelBaseType::ComplexOrRealType ComplexOrRealType;
	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename PsimagLite::Matrix<ComplexOrRealType> MatrixType;

	enum {TERM_HOPPING=0,TERM_NINJ=1,TERM_SUPER=2};

	enum {SPIN_UP, SPIN_DOWN};

	static const int FERMION_SIGN = -1;

	HubbardHelper(const GeometryType& geometry,
	              const ModelParamsType& mp,
	              const MatrixType& hoppings,
	              const bool& hasJcoupling,
	              const bool& hasCoulombCoupling)
	    : geometry_(geometry),
	      mp_(mp),
	      hoppings_(hoppings),
	      hasJcoupling_(hasJcoupling),
	      hasCoulombCoupling_(hasCoulombCoupling)
	{}

	//! Gf. related functions below:
	void setupHamiltonian(SparseMatrixType& matrix,
	                      const BasisBaseType& basis) const
	{
		SizeType hilbert=basis.size();
		typename PsimagLite::Vector<RealType>::Type diag(hilbert);
		calcDiagonalElements(diag, basis);

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

	void matrixVectorProduct(VectorType &x,
	                         VectorType const &y,
	                         const BasisBaseType& basis) const
	{
		SizeType hilbert=basis.size();
		typename PsimagLite::Vector<RealType>::Type diag(hilbert);
		calcDiagonalElements(diag,basis);
		for (SizeType ispace=0;ispace<hilbert;ispace++)
			x[ispace] += diag[ispace]*y[ispace];
		diag.clear();

		SizeType nsite = geometry_.numberOfSites();
		// Calculate off-diagonal elements AND store matrix
		for (SizeType ispace=0;ispace<hilbert;ispace++) {
			SparseRowType sparseRow;
			WordType ket1 = basis(ispace,SPIN_UP);
			WordType ket2 = basis(ispace,SPIN_DOWN);
			for (SizeType i=0;i<nsite;i++) {
				setHoppingTerm(sparseRow,ket1,ket2,i,basis);
				setJTermOffDiagonal(sparseRow,ket1,ket2,i,basis);
			}
			x[ispace] += sparseRow.finalize(y);
		}
	}

private:

	void calcDiagonalElements(typename PsimagLite::Vector<RealType>::Type& diag,
	                          const BasisBaseType& basis) const
	{
		const RealType zeroPointFive = 0.5;
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
					if (PsimagLite::real(value) == 0 && PsimagLite::imag(value) == 0) continue;
					s += value*zeroPointFive* // double counting i,j
					        szTerm(ket1,ket2,i,basis)*
					        szTerm(ket1,ket2,j,basis);
				}

				// Coulomb
				RealType ne = (basis.getN(ket1,ket2,i,SPIN_UP,orb) +
				               basis.getN(ket1,ket2,i,SPIN_DOWN,orb));

				for (SizeType j=0;j<nsite;j++) {
					ComplexOrRealType value = coulombCoupling(i,j);
					if (PsimagLite::real(value) == 0 && PsimagLite::imag(value) == 0) continue;
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

			assert(fabs(PsimagLite::imag(s))<1e-12);
			diag[ispace] = PsimagLite::real(s);
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
			if (PsimagLite::real(h) == 0 && PsimagLite::imag(h) == 0) continue;
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
				if (s1i == 0) cTemp = PsimagLite::conj(cTemp);
				assert(temp<basis.size());
				sparseRow.add(temp,cTemp);
			}

			if (s2i+s2j==1) {
				WordType bra2= ket2 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
				SizeType temp = basis.perfectIndex(ket1,bra2);
				RealType extraSign = (s2i==1) ? FERMION_SIGN : 1;
				RealType tmp2 = basis.doSign(ket1,ket2,i,orb,j,orb,SPIN_DOWN);
				ComplexOrRealType cTemp = h*extraSign*tmp2;
				if (s2j == 0) cTemp = PsimagLite::conj(cTemp);
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
		const RealType zeroPointFive = 0.5;
		for (SizeType j=0;j<geometry_.numberOfSites();j++) {
			ComplexOrRealType value = jCoupling(i,j)*zeroPointFive;
			if (PsimagLite::real(value) == 0 && PsimagLite::imag(value) == 0) continue;
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

	HubbardHelper(const HubbardHelper&) = delete;

	HubbardHelper& operator=(const HubbardHelper&) = delete;

private:

	const GeometryType& geometry_;
	const ModelParamsType& mp_;
	const MatrixType& hoppings_;
	const bool& hasJcoupling_;
    const bool& hasCoulombCoupling_;
};

}
#endif // HUBBARDHELPER_H
