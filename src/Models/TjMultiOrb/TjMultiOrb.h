/*
*/

#ifndef LANCZOS_TJ_1ORB_H
#define LANCZOS_TJ_1ORB_H

#include "CrsMatrix.h"
#include "BasisTjMultiOrbLanczos.h"
#include "BitManip.h"
#include "TypeToString.h"
#include "SparseRow.h"
#include "ParametersTjMultiOrb.h"
#include "ModelBase.h"

namespace LanczosPlusPlus {

template<typename ComplexOrRealType,typename GeometryType,typename InputType>
class Tj1Orb  : public ModelBase<ComplexOrRealType,GeometryType,InputType> {

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef ModelBase<ComplexOrRealType,GeometryType,InputType> BaseType;

	enum {SPIN_UP = ProgramGlobals::SPIN_UP, SPIN_DOWN = ProgramGlobals::SPIN_DOWN};

public:

	typedef ParametersTj1Orb<RealType,InputType> ParametersModelType;
	typedef BasisTj1OrbLanczos<GeometryType> BasisType;
	typedef typename BasisType::BaseType BasisBaseType;
	typedef typename BasisType::WordType WordType;
	typedef typename BaseType::SparseMatrixType SparseMatrixType;
	typedef typename BaseType::VectorType VectorType;
	typedef typename BaseType::VectorSizeType VectorSizeType;
	typedef PsimagLite::SparseRow<SparseMatrixType> SparseRowType;

	static int const FERMION_SIGN = BasisType::FERMION_SIGN;

	Tj1Orb(SizeType nup,
	       SizeType ndown,
	       InputType& io,
	       const GeometryType& geometry)
	    : mp_(io),
	      geometry_(geometry),
	      basis_(geometry,nup,ndown,mp_.orbitals),
	      hoppings_(geometry_.numberOfSites()*geometry_.numberOfSites(),mp_.orbitals*mp_.orbitals),
	      jpm_(geometry_.numberOfSites(),geometry_.numberOfSites()),
	      jzz_(geometry_.numberOfSites(),geometry_.numberOfSites()),
	      w_(geometry_.numberOfSites(),geometry_.numberOfSites())
	{
		SizeType n = geometry_.numberOfSites();

		if (geometry_.terms() != 4) {
			PsimagLite::String msg("Tj1Orb: must have 4 terms\n");
			throw PsimagLite::RuntimeError(msg);
		}

		for (SizeType i=0;i<n;i++) {
			for (SizeType j=0;j<n;j++) {
				for (SizeType orb1 = 0; orb1 < mp_.orbitals; ++orb1) {
					for (SizeType orb2 = 0; orb2 < mp_.orbitals; ++orb2) {
						hoppings_(i+j*n,orb1+orb2*mp_.orbitals) = geometry_(i,orb1,j,orb2,0);
					}
				}

				jpm_(i,j) = geometry_(i,0,j,0,1);
				jzz_(i,j) = geometry_(i,0,j,0,2);
				w_(i,j) = geometry_(i,0,j,0,3);
			}
		}
	}

	~Tj1Orb()
	{
		BaseType::deleteGarbage(garbage_);
	}

	SizeType size() const { return basis_.size(); }

	SizeType orbitals(SizeType) const
	{
		return mp_.orbitals;
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
		typename PsimagLite::Vector<RealType>::Type diag(hilbert,0.0);
		calcDiagonalElements(diag,basis);

		SizeType nsite = geometry_.numberOfSites();

		matrix.resize(hilbert,hilbert);
		// Calculate off-diagonal elements AND store matrix
		SizeType nCounter=0;
		for (SizeType ispace=0;ispace<hilbert;ispace++) {
			SparseRowType sparseRow;
			matrix.setRow(ispace,nCounter);
			WordType ket1 = basis(ispace,SPIN_UP);
			WordType ket2 = basis( ispace,SPIN_DOWN);
			//				std::cout<<"ket1="<<ket1<<" ket2="<<ket2<<"\n";
			// Save diagonal
			sparseRow.add(ispace,diag[ispace]);
			for (SizeType i=0;i<nsite;i++) {
				for (SizeType orb = 0; orb < mp_.orbitals; ++orb) {
					setHoppingTerm(sparseRow,ket1,ket2,i,orb,basis);
					setSplusSminus(sparseRow,ket1,ket2,i,orb,basis);
				}
			}

			nCounter += sparseRow.finalize(matrix);
		}
		matrix.setRow(hilbert,nCounter);
		matrix.checkValidity();
		assert(isHermitian(matrix));
	}

	bool hasNewParts(std::pair<SizeType,SizeType>& newParts,
	                 SizeType what,
	                 SizeType spin,
	                 const PairType& orbs) const
	{
		if (what==ProgramGlobals::OPERATOR_C || what==ProgramGlobals::OPERATOR_CDAGGER)
			return hasNewPartsCorCdagger(newParts,what,spin,orbs);

		if (what == ProgramGlobals::OPERATOR_SPLUS ||
		        what == ProgramGlobals::OPERATOR_SMINUS)
			return hasNewPartsSplusOrMinus(newParts,what,spin,orbs);

		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) +  "\n";
		str += PsimagLite::String("hasNewParts: unsupported operator ");
		str += ProgramGlobals::id2Operator(what) + "\n";
		throw std::runtime_error(str.c_str());
	}

	const GeometryType& geometry() const { return geometry_; }

	const BasisType& basis() const { return basis_; }

	void printBasis(std::ostream& os) const
	{
		os<<basis_;
	}

	PsimagLite::String name() const { return __FILE__; }

	BasisType* createBasis(SizeType nup, SizeType ndown) const
	{
		BasisType* ptr = new BasisType(geometry_,nup,ndown,mp_.orbitals);
		garbage_.push_back(ptr);
		return ptr;
	}

	void print(std::ostream& os) const { os<<mp_; }

	void printOperators(std::ostream& os) const
	{
		SizeType nup = basis_.electrons(SPIN_UP);
		SizeType ndown = basis_.electrons(SPIN_DOWN);
		os<<"#SectorSource 2 "<<nup<<" "<<ndown<<"\n";
		SizeType spin = SPIN_UP;
		for (SizeType site = 0; site < geometry_.numberOfSites(); ++site)
			printOperatorC(site,spin,os);
	}

private:

	void printOperatorC(SizeType site, SizeType spin, std::ostream& os) const
	{
		SizeType nup = basis_.electrons(SPIN_UP);
		SizeType ndown = basis_.electrons(SPIN_DOWN);
		if (nup == 0) {
			os<<"#Operator_c_"<<spin<<"_"<<site<<"\n";
			os<<"#SectorDest 0\n"; //bogus
			os<<"#Matrix\n";
			os<<"0 0\n";
			return;
		}

		BasisType*  basis = createBasis(nup-1, ndown);
		VectorSizeType opt(2,0);
		opt[0] = site;
		opt[1] = spin;
		MatrixType matrix;
		setupOperator(matrix,*basis,"c",opt);
		os<<"#Operator_c_"<<spin<<"_"<<site<<"\n";
		os<<"#SectorDest 2 "<<(nup-1)<<" "<<ndown<<"\n";
		os<<"#Matrix\n";
		os<<matrix;
	}

	void setupOperator(MatrixType& matrix,
	                   const BasisBaseType& basis,
	                   PsimagLite::String operatorName,
	                   const VectorSizeType& operatorOptions) const
	{
		SizeType hilbertDest = basis.size();
		SizeType hilbertSrc = basis_.size();
		SizeType nsite = geometry_.numberOfSites();
		SizeType id = 0;
		if (operatorName == "c") {
			id = ProgramGlobals::OPERATOR_C;
		} else {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "operator " + operatorName + " is unimplemented for this model\n";
			throw PsimagLite::RuntimeError(str);
		}

		if (operatorOptions.size() < 2) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "operator "+ operatorName + " needs at least two options\n";
			throw PsimagLite::RuntimeError(str);
		}

		SizeType site = operatorOptions[0];
		if (site >= nsite) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "site requested " + ttos(site);
			str += " but number of sites= " + ttos(nsite) + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		SizeType spin = operatorOptions[1];
		matrix.resize(hilbertSrc,hilbertDest);
		matrix.setTo(0.0);
		SizeType orb = 0;

		for (SizeType ispace=0;ispace<hilbertSrc;ispace++) {
			WordType ket1 = basis_(ispace,SPIN_UP);
			WordType ket2 = basis_(ispace,SPIN_DOWN);
			WordType bra = ket1;
			// assumes OPERATOR_C
			bool b = basis.getBra(bra,ket1,ket2,id,site,spin);
			if (!b) continue;
			SizeType index = basis.perfectIndex(bra,ket2);

			matrix(ispace,index) = basis.doSignGf(bra,ket2,site,spin,orb);
		}
	}

	bool hasNewPartsCorCdagger(std::pair<SizeType,SizeType>& newParts,
	                           SizeType what,
	                           SizeType spin,
	                           const PairType&) const
	{
		int newPart1=basis_.electrons(SPIN_UP);
		int newPart2=basis_.electrons(SPIN_DOWN);
		int c = (what==ProgramGlobals::OPERATOR_CDAGGER) ? 1 : -1;
		if (spin==SPIN_UP) newPart1 += c;
		else newPart2 += c;

		if (newPart1<0 || newPart2<0) return false;
		SizeType nsite = geometry_.numberOfSites();
		if (SizeType(newPart1)>nsite || SizeType(newPart2)>nsite) return false;
		if (newPart1==0 && newPart2==0) return false;
		if (SizeType(newPart1+newPart2)>nsite) return false; // no double occupancy
		newParts.first = SizeType(newPart1);
		newParts.second = SizeType(newPart2);
		return true;
	}

	bool hasNewPartsSplusOrMinus(std::pair<SizeType,SizeType>& newParts,
	                             SizeType what,
	                             SizeType spin,
	                             const PairType&) const
	{
		int newPart1=basis_.electrons(SPIN_UP);
		int newPart2=basis_.electrons(SPIN_DOWN);
		int c = (what==ProgramGlobals::OPERATOR_SPLUS) ? 1 : -1;
		if (spin==SPIN_UP) {
			newPart1 += c;
			newPart2 -= c;
		} else {
			newPart2 += c;
			newPart1 -= c;
		}

		if (newPart1<0 || newPart2<0) return false;

		SizeType nsite = geometry_.numberOfSites();
		if (SizeType(newPart1)>nsite || SizeType(newPart2)>nsite) return false;
		if (newPart1==0 && newPart2==0) return false;
		if (SizeType(newPart1+newPart2)>nsite) return false; // no double occupancy
		newParts.first = SizeType(newPart1);
		newParts.second = SizeType(newPart2);
		return true;
	}

	void calcDiagonalElements(typename PsimagLite::Vector<RealType>::Type& diag,
	                          const BasisBaseType& basis) const
	{
		SizeType hilbert=basis.size();
		SizeType nsite = geometry_.numberOfSites();

		// Calculate diagonal elements
		for (SizeType ispace=0;ispace<hilbert;ispace++) {
			WordType ket1 = basis(ispace,SPIN_UP);
			WordType ket2 = basis(ispace,SPIN_DOWN);
			ComplexOrRealType s=0;
			for (SizeType i=0;i<nsite;i++) {
				for (SizeType orb = 0; orb < mp_.orbitals; ++orb) {
					int niup = basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP,orb);
					int nidown = basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN,orb);

					if (i < mp_.potentialV.size()) {
						s += mp_.potentialV[i+orb*nsite]*niup;
						s += mp_.potentialV[i+orb*nsite+mp_.orbitals*nsite]*nidown;
					}

					for (SizeType j=i+1;j<nsite;j++) {
						for (SizeType orb2 = 0; orb2 < mp_.orbitals; ++orb2) {
							int njup = basis.isThereAnElectronAt(ket1,ket2,j,SPIN_UP,orb2);
							int njdown = basis.isThereAnElectronAt(ket1,ket2,j,SPIN_DOWN,orb2);

							// Sz Sz term:
							s += (niup-nidown) * (njup - njdown)  * jzz_(i,j)*0.25;

							// ni nj term
							s+= (niup+nidown) * (njup + njdown) * w_(i,j);
						}
					}
				}
			}

			assert(fabs(std::imag(s))<1e-12);
			diag[ispace] = std::real(s);
		}
	}

	void setHoppingTerm(SparseRowType &sparseRow,
	                    const WordType& ket1,
	                    const WordType& ket2,
	                    SizeType i,
	                    SizeType orb,
	                    const BasisBaseType &basis) const
	{
		WordType s1i=(ket1 & BasisType::bitmask(i*mp_.orbitals+orb));
		if (s1i>0) s1i=1;
		WordType s2i=(ket2 & BasisType::bitmask(i*mp_.orbitals+orb));
		if (s2i>0) s2i=1;

		SizeType nsite = geometry_.numberOfSites();

		// Hopping term
		for (SizeType j=0;j<nsite;j++) {
			for (SizeType orb2 = 0; orb2 < mp_.orbitals; ++orb2) {
				if (j<i) continue;
				ComplexOrRealType h = hoppings_(i+j*nsite,orb+orb2*mp_.orbitals);
				if (std::real(h) == 0 && std::imag(h) == 0) continue;
				WordType s1j= (ket1 & BasisType::bitmask(j*mp_.orbitals+orb2));
				if (s1j>0) s1j=1;
				WordType s2j= (ket2 & BasisType::bitmask(j*mp_.orbitals+orb2));
				if (s2j>0) s2j=1;

				if (s1i+s1j==1 && !(s1j==0 && s2j>0) && !(s1j>0 && s2i>0)) {
					WordType bra1= ket1 ^ (BasisType::bitmask(i*mp_.orbitals+orb)|
					                       BasisType::bitmask(j*mp_.orbitals+orb2));
					SizeType temp = basis.perfectIndex(bra1,ket2);
					RealType extraSign = (s1i==1) ? FERMION_SIGN : 1;
					RealType tmp2 = basis_.doSign(ket1,ket2,i,orb,j,orb2,SPIN_UP);
					ComplexOrRealType cTemp = h*extraSign*tmp2;
					sparseRow.add(temp,cTemp);
				}

				if (s2i+s2j==1 && !(s2j==0 && s1j>0) && !(s2j>0 && s1i>0)) {
					WordType bra2= ket2 ^(BasisType::bitmask(i*mp_.orbitals+orb)|
					                      BasisType::bitmask(j*mp_.orbitals+orb2));
					SizeType temp = basis.perfectIndex(ket1,bra2);
					RealType extraSign = (s2i==1) ? FERMION_SIGN : 1;
					RealType tmp2 = basis_.doSign(ket1,ket2,i,orb,j,orb2,SPIN_DOWN);
					ComplexOrRealType cTemp = h*extraSign*tmp2;
					sparseRow.add(temp,cTemp);
				}
			}
		}
	}

	void setSplusSminus(SparseRowType &sparseRow,
	                    const WordType& ket1,
	                    const WordType& ket2,
	                    SizeType i,
	                    SizeType orb,
	                    const BasisBaseType &basis) const
	{
		WordType s1i=(ket1 & BasisType::bitmask(i*mp_.orbitals+orb));
		if (s1i>0) s1i=1;
		WordType s2i=(ket2 & BasisType::bitmask(i*mp_.orbitals+orb));
		if (s2i>0) s2i=1;

		SizeType nsite = geometry_.numberOfSites();

		// Hopping term
		for (SizeType j=0;j<nsite;j++) {
			for (SizeType orb2 = 0; orb2 < mp_.orbitals; ++orb2) {
				if (j<i) continue;
				ComplexOrRealType h = jpm_(i,j)*0.5;
				if (std::real(h) == 0 && std::imag(h) == 0) continue;
				WordType s1j= (ket1 & BasisType::bitmask(j*mp_.orbitals+orb2));
				if (s1j>0) s1j=1;
				WordType s2j= (ket2 & BasisType::bitmask(j*mp_.orbitals+orb2));
				if (s2j>0) s2j=1;

				if (s1i==1 && s1j==0 && s2i==0 && s2j==1) {
					WordType bra1= ket1 ^ BasisType::bitmask(i*mp_.orbitals+orb);
					bra1 |= BasisType::bitmask(j*mp_.orbitals+orb2);
					WordType bra2= ket2 | BasisType::bitmask(i*mp_.orbitals+orb);
					bra2 ^= BasisType::bitmask(j*mp_.orbitals+orb2);
					SizeType temp = basis.perfectIndex(bra1,bra2);
					sparseRow.add(temp,h*signSplusSminus(i*mp_.orbitals+orb,
					                                     j*mp_.orbitals+orb2,
					                                     bra1,
					                                     bra2));
				}

				if (s1i==0 && s1j==1 && s2i==1 && s2j==0) {
					WordType bra1= ket1 | BasisType::bitmask(i*mp_.orbitals+orb);
					bra1 ^= BasisType::bitmask(j*mp_.orbitals+orb2);
					WordType bra2= ket2 ^ BasisType::bitmask(i*mp_.orbitals+orb);
					bra2 |= BasisType::bitmask(j*mp_.orbitals+orb2);
					SizeType temp = basis.perfectIndex(bra1,bra2);
					sparseRow.add(temp,h*signSplusSminus(i*mp_.orbitals+orb,
					                                     j*mp_.orbitals+orb2,
					                                     bra1,
					                                     bra2));
				}
			}
		}
	}

	RealType signSplusSminus(SizeType i,
	                         SizeType j,
	                         const WordType& bra1,
	                         const WordType& bra2) const
	{
		int s = 1;
		if (j > 0) s *= parityFrom(0,j-1,bra2);
		if (i > 0) s *= parityFrom(0,i-1,bra2);
		if (i > 0) s *= parityFrom(0,i-1,bra1);
		if (j > 0) s *= parityFrom(0,j-1,bra1);
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
	PsimagLite::Matrix<ComplexOrRealType> hoppings_;
	PsimagLite::Matrix<ComplexOrRealType> jpm_;
	PsimagLite::Matrix<ComplexOrRealType> jzz_;
	PsimagLite::Matrix<ComplexOrRealType> w_;
	mutable typename PsimagLite::Vector<BasisType*>::Type garbage_;
}; // class Tj1Orb
} // namespace LanczosPlusPlus
#endif

