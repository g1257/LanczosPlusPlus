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
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelBaseType::RahulOperatorType RahulOperatorType;
	typedef typename ModelBaseType::VectorRahulOperatorType VectorRahulOperatorType;

	enum {SPIN_UP = ProgramGlobals::SPIN_UP, SPIN_DOWN = ProgramGlobals::SPIN_DOWN};

	enum TermEnum {HOPPING = 0, NINJ = 1, SUPER = 2};

	static const int FERMION_SIGN = -1;

	HubbardHelper(const GeometryType& geometry,
	              const ModelParamsType& mp)
	    : geometry_(geometry),
	      mp_(mp),
	      hoppings_(geometry.numberOfSites(), geometry.numberOfSites()),
	      hasJcoupling_(mp_.model == "SuperHubbardExtended"),
	      hasCoulombCoupling_(mp_.model == "HubbardOneBandExtended" ||
	                          mp_.model == "SuperHubbardExtended"),
	      hasRashba_(mp_.model == "HubbardOneBandRashbaSOC")
	{
		const bool hasSpinOrbitKaneMele = (mp_.model == "KaneMeleHubbard");

		if (hasCoulombCoupling_ && geometry_.terms() < 2)
			err("HubbardHelper::ctor(): ColoumbCoupling\n");

		if (hasJcoupling_ && geometry_.terms()<3)
			err("HubbardHelper::ctor(): jCoupling\n");

		if (hasSpinOrbitKaneMele && geometry_.terms() != 2)
			err("HubbardHelper::ctor(): KaneMeleHubbard\n");

		if (hasRashba_ && geometry_.terms() != 2)
			err("HubbardHelper::ctor(): Rashba needs two Hamiltonian terms\n");

		const SizeType n = geometry_.numberOfSites();
		if (hasRashba_) rashbaHoppings_.resize(n, n);
		for (SizeType j=0;j<n;j++) {
			for (SizeType i=0;i<j;i++) {

				hoppings_(i, j) = geometry_(i, 0, j, 0, 0);

				if (hasSpinOrbitKaneMele)
					hoppings_(i, j) += geometry_(i, 0, j, 0, 1);

				if (hasRashba_)
					rashbaHoppings_(i, j) = geometry_(i, 0, j, 0, 1);
			}
		}
	}

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

	void rahulMethod(VectorType& psiNew,
	                 const VectorRahulOperatorType& vops,
	                 const VectorSizeType& vsites,
	                 const VectorType& psi,
	                 const BasisBaseType& basis) const
	{
		std::fill(psiNew.begin(), psiNew.end(), 0.0);
		const SizeType hilbert = basis.size();
		const SizeType nops = vops.size();
		for (SizeType ispace = 0; ispace < hilbert; ++ispace) {
			WordType ket1 = basis(ispace, SPIN_UP);
			WordType ket2 = basis(ispace, SPIN_DOWN);
			assert(ispace < psi.size());
			ComplexOrRealType value = psi[ispace];
			WordType ketp1 = ket1;
			WordType ketp2 = ket2;
			bool nonZero = false;
			for (SizeType jj = 0; jj < nops; ++jj) {
				const SizeType j = nops - jj - 1; // start from the end
				const RahulOperatorType& op = vops[j];
				assert(j < vsites.size());
				const SizeType site = vsites[j];

				ComplexOrRealType result = 0;
				nonZero = (op.dof() == SPIN_UP)
				        ? applyOperator(ketp1, result, op, site)
				        : applyOperator(ketp2, result, op, site);
				if (!nonZero) break;

				if (op.isFermionic()) {
					const SizeType overUp = (op.dof() == SPIN_UP)
					        ? 0
					        : PsimagLite::BitManip::count(ket1);

					if (overUp & 1) result *= FERMION_SIGN;

					const WordType ket1or2 = (op.dof() == SPIN_UP) ? ketp1 : ketp2;
					const ComplexOrRealType fsign = ProgramGlobals::doSign(ket1or2, site);
					result *= fsign;
				}

				value *= result;
			}

			if (!nonZero) continue;
			SizeType newI = basis.perfectIndex(ketp1, ketp2);
			assert(newI < psiNew.size());
			psiNew[newI] += value;
		}
	}

private:

	bool applyOperator(WordType& ketp,
	                   ComplexOrRealType& result,
	                   const RahulOperatorType& op,
	                   SizeType site) const
	{
		const WordType ket = ketp;
		const WordType mask = BasisType::bitmask(site);
		WordType s = (ket & mask);
		bool sbit = (s > 0);
		bool sbitSaved = sbit;
		bool nonZero = op.actOn(sbit, result);
		if (!nonZero) return false;
		if (sbitSaved != sbit) ketp ^= mask;
		return true;
	}

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
					ComplexOrRealType value = 0.5*coulombCoupling(i,j);
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

		const SizeType nsite = geometry_.numberOfSites();
		const SizeType orb = 0;

		// Hopping term
		for (SizeType j=0;j<nsite;j++) {
			if (j < i) continue;
			const ComplexOrRealType& h = hoppings_(i,j);
			const bool hasHop = (PsimagLite::real(h) != 0 || PsimagLite::imag(h) != 0);
			WordType s1j= (ket1 & BasisType::bitmask(j));
			if (s1j>0) s1j=1;
			WordType s2j= (ket2 & BasisType::bitmask(j));
			if (s2j>0) s2j=1;

			if (hasHop && s1i + s1j == 1) {
				WordType bra1= ket1 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
				SizeType temp = basis.perfectIndex(bra1,ket2);
				RealType extraSign = (s1i==1) ? FERMION_SIGN : 1;
				RealType tmp2 = basis.doSign(ket1,ket2,i,orb,j,orb,SPIN_UP);
				ComplexOrRealType cTemp = h*extraSign*tmp2;
				if (s1i == 0) cTemp = PsimagLite::conj(cTemp);
				assert(temp<basis.size());
				sparseRow.add(temp,cTemp);
			}

			if (hasHop && s2i + s2j == 1) {
				WordType bra2= ket2 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
				SizeType temp = basis.perfectIndex(ket1,bra2);
				RealType extraSign = (s2i==1) ? FERMION_SIGN : 1;
				RealType tmp2 = basis.doSign(ket1,ket2,i,orb,j,orb,SPIN_DOWN);
				ComplexOrRealType cTemp = h*extraSign*tmp2;
				if (s2j == 0) cTemp = PsimagLite::conj(cTemp);
				assert(temp<basis.size());
				sparseRow.add(temp,cTemp);
			}

			if (!hasRashba_) continue;

			const ComplexOrRealType& hr = rashbaHoppings_(i, j);
			if (PsimagLite::real(hr) == 0 && PsimagLite::imag(hr) == 0) continue;

			if (s1i + s2j == 1) {
				WordType bra1= ket1 ^ (BasisType::bitmask(i));
				WordType bra2= ket2 ^ (BasisType::bitmask(j));
				SizeType temp = basis.perfectIndex(bra1, bra2);
				RealType extraSign = (s2j == 0) ? FERMION_SIGN : 1;
				RealType tmp2 = ProgramGlobals::doSign(ket1, i)*ProgramGlobals::doSign(ket2, j);
				const SizeType count1 = PsimagLite::BitManip::count(ket1); // + s1i;
				if (count1 & 1) tmp2 *= FERMION_SIGN;
				ComplexOrRealType cTemp = hr*extraSign*tmp2;
				if (s2j == 0) cTemp = PsimagLite::conj(cTemp);
				assert(temp<basis.size());
				sparseRow.add(temp,cTemp);
			}

			if (s2i + s1j == 1) {
				WordType bra1= ket1 ^ (BasisType::bitmask(j));
				WordType bra2= ket2 ^(BasisType::bitmask(i));
				SizeType temp = basis.perfectIndex(bra1, bra2);
				RealType extraSign = (s2i == 0) ? FERMION_SIGN : 1;
				RealType tmp2 = ProgramGlobals::doSign(ket1, j)*ProgramGlobals::doSign(ket2, i);
				const SizeType count1 = PsimagLite::BitManip::count(ket1); // + s1j;
				if (count1 & 1) tmp2 *= FERMION_SIGN;
				ComplexOrRealType cTemp = hr*extraSign*tmp2;
				if (s2i == 0) cTemp = PsimagLite::conj(cTemp);
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
		return geometry_(i, 0, j, 0, TermEnum::SUPER);
	}

	ComplexOrRealType coulombCoupling(SizeType i,SizeType j) const
	{
		if (!hasCoulombCoupling_) return 0;
		return geometry_(i, 0, j, 0, TermEnum::NINJ);
	}

	HubbardHelper(const HubbardHelper&) = delete;

	HubbardHelper& operator=(const HubbardHelper&) = delete;

private:

	const GeometryType& geometry_;
	const ModelParamsType& mp_;
	MatrixType hoppings_;
	MatrixType rashbaHoppings_;
	bool hasJcoupling_;
	bool hasCoulombCoupling_;
	bool hasRashba_;
};

}
#endif // HUBBARDHELPER_H
