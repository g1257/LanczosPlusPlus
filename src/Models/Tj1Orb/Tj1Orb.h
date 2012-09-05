
/*
*/

#ifndef TJ_1ORB_H
#define TJ_1ORB_H

#include "CrsMatrix.h"
#include "BasisTj1OrbLanczos.h"
#include "BitManip.h"
#include "TypeToString.h"
#include "SparseRow.h"

namespace LanczosPlusPlus {

	template<typename RealType_,typename ParametersType,typename GeometryType_>
	class Tj1Orb {

		typedef PsimagLite::Matrix<RealType_> MatrixType;

	public:

		typedef GeometryType_ GeometryType;
		typedef PsimagLite::CrsMatrix<RealType_> SparseMatrixType;
		typedef PsimagLite::SparseRow<SparseMatrixType> SparseRowType;
		typedef BasisTj1OrbLanczos<GeometryType> BasisType;
		typedef typename BasisType::WordType WordType;
		typedef RealType_ RealType;
		typedef std::vector<RealType> VectorType;

		enum {SPIN_UP=BasisType::SPIN_UP,SPIN_DOWN=BasisType::SPIN_DOWN};

		enum {DESTRUCTOR=BasisType::DESTRUCTOR,CONSTRUCTOR=BasisType::CONSTRUCTOR};

		static int const FERMION_SIGN = BasisType::FERMION_SIGN;

		Tj1Orb(size_t nup,
		                  size_t ndown,
		                  const ParametersType& mp,
		                  GeometryType& geometry)
		: mp_(mp),
		  geometry_(geometry),
		  basis_(geometry,nup,ndown),
		  hoppings_(geometry_.numberOfSites(),geometry_.numberOfSites()),
		  j_(geometry_.numberOfSites(),geometry_.numberOfSites()),
		  w_(geometry_.numberOfSites(),geometry_.numberOfSites())
		{
			size_t n = geometry_.numberOfSites();
			for (size_t i=0;i<n;i++) {
				for (size_t j=0;j<n;j++) {
					hoppings_(i,j) = geometry_(i,0,j,0,0);
					j_(i,j) = geometry_(i,0,j,0,1);
					w_(i,j) = geometry_(i,0,j,0,2);
				}
			}
		}

		size_t size() const { return basis_.size(); }

		size_t orbitals() const { return 1; }

		void setupHamiltonian(SparseMatrixType &matrix) const
		{
			setupHamiltonian(matrix,basis_);
		}

		//! Gf. related functions below:
		void setupHamiltonian(SparseMatrixType &matrix,
		                      const BasisType &basis) const
		{
			size_t hilbert=basis.size();
			std::vector<RealType> diag(hilbert,0.0);
			calcDiagonalElements(diag,basis);
			
			size_t nsite = geometry_.numberOfSites();

			matrix.resize(hilbert,hilbert);
			// Calculate off-diagonal elements AND store matrix
			size_t nCounter=0;
			for (size_t ispace=0;ispace<hilbert;ispace++) {
				SparseRowType sparseRow;
				matrix.setRow(ispace,nCounter);
				WordType ket1 = basis(ispace,SPIN_UP);
				WordType ket2 = basis(ispace,SPIN_DOWN);
//				std::cout<<"ket1="<<ket1<<" ket2="<<ket2<<"\n";
				// Save diagonal
				sparseRow.add(ispace,diag[ispace]);
				for (size_t i=0;i<nsite;i++) {
					setHoppingTerm(sparseRow,ket1,ket2,i,basis);
					setSplusSminus(sparseRow,ket1,ket2,i,basis);
				}
				nCounter += sparseRow.finalize(matrix);
			}
			matrix.setRow(hilbert,nCounter);
			matrix.checkValidity();
			assert(isHermitian(matrix));
//			PsimagLite::Matrix<RealType> m;
//			crsMatrixToFullMatrix(m,matrix);
//			std::cout<<m;
		}

		bool hasNewParts(std::pair<size_t,size_t>& newParts,
		                 size_t type,
				 size_t spin,
				 const std::pair<size_t,size_t>& orbs) const
		{
			int newPart1=basis_.electrons(SPIN_UP);
			int newPart2=basis_.electrons(SPIN_DOWN);
			int c = (type&1) ? 1 : -1;
			if (spin==SPIN_UP) newPart1 += c;
			else newPart2 += c;

			if (newPart1<0 || newPart2<0) return false;
			size_t nsite = geometry_.numberOfSites();
			if (size_t(newPart1)>nsite || size_t(newPart2)>nsite) return false;
			if (newPart1==0 && newPart2==0) return false;
			if (size_t(newPart1+newPart2)>nsite) return false; // no double occupancy
			newParts.first = size_t(newPart1);
			newParts.second = size_t(newPart2);
			return true;
		}

		void getModifiedState(std::vector<RealType>& modifVector,
		                      const std::vector<RealType>& gsVector,
		                      const BasisType& basisNew,
		                      size_t type,
		                      size_t isite,
		                      size_t jsite,
		                      size_t spin) const
		{
			size_t what= (type&1) ? CONSTRUCTOR : DESTRUCTOR;

			modifVector.resize(basisNew.size());
			for (size_t temp=0;temp<modifVector.size();temp++)
				modifVector[temp]=0.0;

			size_t orb = 0; // bogus orbital index, no orbitals in this model
			accModifiedState(modifVector,basisNew,gsVector,what,isite,spin,orb,1);
			std::cerr<<"isite="<<isite<<" type="<<type;
			std::cerr<<" modif="<<(modifVector*modifVector)<<"\n";
			if (isite==jsite) return;

			int isign= (type>1) ? -1 : 1;
			accModifiedState(modifVector,basisNew,gsVector,what,jsite,spin,orb,isign);
			std::cerr<<"jsite="<<jsite<<" type="<<type;
			std::cerr<<" modif="<<(modifVector*modifVector)<<"\n";
		}

		const GeometryType& geometry() const { return geometry_; }

		const BasisType& basis() const { return basis_; }

		//! Gf Related functions:
		void accModifiedState(std::vector<RealType> &z,
		const BasisType& newBasis,
		const std::vector<RealType> &gsVector,
		size_t what,
		size_t site,
		size_t spin,
		size_t orb,
		int isign) const
		{
			for (size_t ispace=0;ispace<basis_.size();ispace++) {
				WordType ket1 = basis_(ispace,SPIN_UP);
				WordType ket2 = basis_(ispace,SPIN_DOWN);
				int temp = newBasis.getBraIndex(ket1,ket2,what,site,spin);
// 				int temp= getBraIndex(mysign,ket1,ket2,newBasis,what,site,spin);
				if (temp>=0 && size_t(temp)>=z.size()) {
					std::string s = "old basis=" + ttos(basis_.size());
					s += " newbasis=" + ttos(newBasis.size());
					s += "\n";
					s += "what=" + ttos(what) + " spin=" + ttos(spin);
					s += " site=" + ttos(site);
					s += "ket1=" + ttos(ket1) + " and ket2=" + ttos(ket2);
					s += "\n";
					s += "getModifiedState: z.size=" + ttos(z.size());
					s += " but temp=" + ttos(temp) + "\n";
					throw std::runtime_error(s.c_str());
				}
				if (temp<0) continue;
				int mysign = basis_.doSignGf(ket1,ket2,site,spin);
				z[temp] += isign*mysign*gsVector[ispace];
			}
		}

		void printBasis(std::ostream& os) const
		{
			os<<basis_;
		}

	private:

		void calcDiagonalElements(std::vector<RealType>& diag,
		                          const BasisType &basis) const
		{
			size_t hilbert=basis.size();
			size_t nsite = geometry_.numberOfSites();

			// Calculate diagonal elements
			for (size_t ispace=0;ispace<hilbert;ispace++) {
				WordType ket1 = basis(ispace,SPIN_UP);
				WordType ket2 = basis(ispace,SPIN_DOWN);
				RealType s=0;
				for (size_t i=0;i<nsite;i++) {
					for (size_t j=i+1;j<nsite;j++) {

						int niup = basis.isThereAnElectronAt(ket1,ket2,i,SPIN_UP);
						int nidown = basis.isThereAnElectronAt(ket1,ket2,i,SPIN_DOWN);
						int njup = basis.isThereAnElectronAt(ket1,ket2,j,SPIN_UP);
						int njdown = basis.isThereAnElectronAt(ket1,ket2,j,SPIN_DOWN);

						// Sz Sz term:
						s += (niup-nidown) * (njup - njdown)  * j_(i,j)*0.25;

						// ni nj term
						s+= (niup+nidown) * (njup + njdown) * w_(i,j);
					}
				}
				diag[ispace]=s;
			}
		}

		void setHoppingTerm(SparseRowType &sparseRow,
		                    const WordType& ket1,
		                    const WordType& ket2,
		                    size_t i,
		                    const BasisType &basis) const
		{
			WordType s1i=(ket1 & BasisType::bitmask(i));
			if (s1i>0) s1i=1;
			WordType s2i=(ket2 & BasisType::bitmask(i));
			if (s2i>0) s2i=1;

			size_t nsite = geometry_.numberOfSites();

			// Hopping term
			for (size_t j=0;j<nsite;j++) {
				if (j<i) continue;
				RealType h = hoppings_(i,j);
				if (h==0) continue;
				WordType s1j= (ket1 & BasisType::bitmask(j));
				if (s1j>0) s1j=1;
				WordType s2j= (ket2 & BasisType::bitmask(j));
				if (s2j>0) s2j=1;

				if (s1i+s1j==1 && !(s1j==0 && s2j>0) && !(s1j>0 && s2i>0)) {
					WordType bra1= ket1 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
					size_t temp = basis.perfectIndex(bra1,ket2);
					int extraSign = (s1i==1) ? FERMION_SIGN : 1;
					RealType cTemp = h*extraSign*basis_.doSign(ket1,ket2,i,j,SPIN_UP);
					sparseRow.add(temp,cTemp);
				}

				if (s2i+s2j==1 && !(s2j==0 && s1j>0) && !(s2j>0 && s1i>0)) {
					WordType bra2= ket2 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
					size_t temp = basis.perfectIndex(ket1,bra2);
					int extraSign = (s2i==1) ? FERMION_SIGN : 1;
					RealType cTemp = h*extraSign*basis_.doSign(ket1,ket2,i,j,SPIN_DOWN);
					sparseRow.add(temp,cTemp);
				}
			}
		}

		void setSplusSminus(SparseRowType &sparseRow,
							const WordType& ket1,
							const WordType& ket2,
							size_t i,
							const BasisType &basis) const
		{
			WordType s1i=(ket1 & BasisType::bitmask(i));
			if (s1i>0) s1i=1;
			WordType s2i=(ket2 & BasisType::bitmask(i));
			if (s2i>0) s2i=1;

			size_t nsite = geometry_.numberOfSites();

			// Hopping term
			for (size_t j=0;j<nsite;j++) {
				if (j<i) continue;
				RealType h = j_(i,j)*0.5;
				if (h==0) continue;
				WordType s1j= (ket1 & BasisType::bitmask(j));
				if (s1j>0) s1j=1;
				WordType s2j= (ket2 & BasisType::bitmask(j));
				if (s2j>0) s2j=1;

				if (s1i==1 && s1j==0 && s2i==0 && s2j==1) {
					WordType bra1= ket1 ^ BasisType::bitmask(i);
					bra1 |= BasisType::bitmask(j);
					WordType bra2= ket2 | BasisType::bitmask(i);
					bra2 ^= BasisType::bitmask(j);
					size_t temp = basis.perfectIndex(bra1,bra2);
					sparseRow.add(temp,h*signSplusSminus(i,j,bra1,bra2));
				}

				if (s1i==0 && s1j==1 && s2i==1 && s2j==0) {
					WordType bra1= ket1 | BasisType::bitmask(i);
					bra1 ^= BasisType::bitmask(j);
					WordType bra2= ket2 ^ BasisType::bitmask(i);
					bra2 |= BasisType::bitmask(j);
					size_t temp = basis.perfectIndex(bra1,bra2);
					sparseRow.add(temp,h*signSplusSminus(i,j,bra1,bra2));
				}
			}
		}

		int signSplusSminus(size_t i, size_t j,const WordType& bra1, const WordType& bra2) const
		{
			size_t n = geometry_.numberOfSites();
			int s = 1.0;
			if (j>0) s *= parityFrom(0,j-1,bra2);
			if (i>0) s *= parityFrom(0,i-1,bra2);
			if (i<n-1) s *= parityFrom(i+1,n-1,bra1);
			if (j<n-1) s *= parityFrom(j+1,n-1,bra1);
			return -s;
		}

		// from i to j including i and j
		// assumes i<=j
		int parityFrom(size_t i,size_t j,const WordType& ket) const
		{
			if (i==j) return (BasisType::bitmask(j) & ket) ? -1 : 1;
			assert(i<j);
			//j>i>=0 now
			WordType mask = ket;
			mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
			int s=(PsimagLite::BitManip::count(mask) & 1) ? -1 : 1; // Parity of up between i and j
			// Is there something of this species at i?
			if (BasisType::bitmask(i) & ket) s = -s;
			// Is there something of this species at j?
			if (BasisType::bitmask(j) & ket) s = -s;
			return s;
		}

		const ParametersType& mp_;
		const GeometryType& geometry_;
		BasisType basis_;
		PsimagLite::Matrix<RealType> hoppings_;
		PsimagLite::Matrix<RealType> j_;
		PsimagLite::Matrix<RealType> w_;

	}; // class Tj1Orb
} // namespace LanczosPlusPlus
#endif

