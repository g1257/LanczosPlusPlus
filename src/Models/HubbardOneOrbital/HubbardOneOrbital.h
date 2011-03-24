
/*
*/

#ifndef HUBBARDLANCZOS_H
#define HUBBARDLANCZOS_H

#include "CrsMatrix.h"
#include "BasisHubbardLanczos.h"
#include "BitManip.h"

namespace LanczosPlusPlus {
	
	template<typename RealType_,typename ParametersType,
		typename GeometryType>
	class HubbardOneOrbital {
		
		typedef PsimagLite::Matrix<RealType_> MatrixType; 
	public:
		
		typedef BasisHubbardLanczos BasisType;
		typedef typename BasisType::WordType WordType;
		typedef RealType_ RealType;
		typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;
		typedef std::vector<RealType> VectorType;
		enum {SPIN_UP,SPIN_DOWN};
		enum {DESTRUCTOR,CONSTRUCTOR};
		
		
		HubbardOneOrbital(size_t nup,size_t ndown,const ParametersType& mp,GeometryType& geometry)
			: mp_(mp),geometry_(geometry),
			  basis1_(geometry.numberOfSites(),nup),basis2_(geometry.numberOfSites(),ndown)
		{
		}
		
		

		size_t size() const { return basis1_.size()*basis2_.size(); }
		

		void setupHamiltonian(SparseMatrixType &matrix) const
		{
			setupHamiltonian(matrix,basis1_,basis2_);
		}

		//! Calc Green function G(isite,jsite)  (still diagonal in spin)
		template<typename ContinuedFractionCollectionType>
		void greenFunction(
				ContinuedFractionCollectionType& cfCollection,
				int isite,
				int jsite,
				int spin,
				const RealType& shift) const
		{
			typedef typename ContinuedFractionCollectionType::
					ContinuedFractionType ContinuedFractionType;
			for (size_t type=0;type<4;type++) {
				if (isite==jsite && type>1) continue;
				std::vector<RealType> modifVector;
				getModifiedStates(modifVector,type,isite,jsite,spin);
				std::pair<int,int> newParts = getNewParts(type,spin);
				if (newParts.first<0) continue;
				ContinuedFractionType cf;
				calcGf(cf,newParts,modifVector,type,spin);
				cfCollection.push(cf);
			}
		}

	private:
		
		bool getBra(WordType& bra, const WordType& ket,size_t what,size_t i) const
		{
			WordType s1i=(ket & BasisType::bitmask(i));
			if (what==DESTRUCTOR) {
				if (s1i>0) {
					bra = (ket ^ BasisType::bitmask(i));
				} else {
					return false; // cannot destroy, there's nothing
				}
			} else {
				if (s1i==0) {
					bra = (ket ^ BasisType::bitmask(i));
				} else {
					return false; // cannot contruct, there's already one
				}
			}
			return true;
		}
		

		RealType hoppings(size_t i,size_t j) const
		{
			return geometry_(i,0,j,0,0);
		}
		

		void setupHamiltonian(SparseMatrixType &matrix,const BasisType &basis1,const BasisType& basis2) const
		{
			// Calculate diagonal elements AND count non-zero matrix elements
			size_t hilbert1=basis1.size();
			size_t hilbert2=basis2.size();
			MatrixType diag(hilbert2,hilbert1);
			size_t nzero = countNonZero(diag,basis1,basis2);
			
			size_t nsite = geometry_.numberOfSites();
			
			// Setup CRS matrix
			matrix.resize(hilbert1*hilbert2,nzero);
			
			// Calculate off-diagonal elements AND store matrix
			size_t nCounter=0;
			matrix.setRow(0,0);
			for (size_t ispace1=0;ispace1<hilbert1;ispace1++) {
				WordType ket1 = basis1[ispace1];
				for (size_t ispace2=0;ispace2<hilbert2;ispace2++) {
					WordType ket2 = basis2[ispace2];
					// Save diagonal
					matrix.setCol(nCounter,ispace2+ispace1*hilbert2);
					RealType cTemp=diag(ispace2,ispace1);
					matrix.setValues(nCounter,cTemp);
					nCounter++;
					for (size_t i=0;i<nsite;i++) {
						WordType s1i=(ket1 & BasisType::bitmask(i));
						WordType s2i=(ket2 & BasisType::bitmask(i));
						if (s1i>0) s1i=1;
						if (s2i>0) s2i=1;
						
						// Hopping term 
						for (size_t j=0;j<nsite;j++) {
							if (j<i) continue;
							RealType tmp = hoppings(i,j);
							if (tmp==0) continue;
							WordType s1j= (ket1 & BasisType::bitmask(j));
							WordType s2j= (ket2 & BasisType::bitmask(j));
							if (s1j>0) s1j=1;
							if (s2j>0) s2j=1;
							if (s1i+s1j==1) {
								WordType bra1= ket1 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
								size_t temp = perfectIndex(basis1,basis2,bra1,ket2);
								matrix.setCol(nCounter,temp);
								cTemp=hoppings(i,j)*doSign(ket1,ket2,i,j,SPIN_UP); // check SIGN FIXME
								if (cTemp==0.0) {
									std::cerr<<"ctemp=0 and hopping="<<hoppings(i,j)<<" and i="<<i<<" and j="<<j<<"\n";
								}
								matrix.setValues(nCounter,cTemp);
								nCounter++;
							}
							if (s2i+s2j==1) {
								WordType bra2= ket2 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
								size_t temp = perfectIndex(basis1,basis2,ket1,bra2);
								matrix.setCol(nCounter,temp);
								cTemp=hoppings(i,j)*doSign(ket1,ket2,i,j,SPIN_DOWN); // Check SIGN FIXME
								matrix.setValues(nCounter,cTemp);
								nCounter++;					
							}
						}
					}
					matrix.setRow(ispace2+hilbert2*ispace1+1,nCounter);
				}
			}
			matrix.setRow(hilbert1*hilbert2,nCounter);
		}
		

		size_t countNonZero(MatrixType& diag,const BasisType &basis1,const BasisType& basis2) const
		{
			size_t hilbert1=basis1.size();
			size_t hilbert2=basis2.size();
			size_t nsite = geometry_.numberOfSites();

			// Calculate diagonal elements AND count non-zero matrix elements
			size_t nzero = 0;
			for (size_t ispace1=0;ispace1<hilbert1;ispace1++) {
				WordType ket1 = basis1[ispace1];
				for (size_t ispace2=0;ispace2<hilbert2;ispace2++) {
					WordType ket2 = basis2[ispace2];
					RealType s=0;
					for (size_t i=0;i<nsite;i++) {
						WordType s1i=(ket1 & BasisType::bitmask(i));
						WordType s2i=(ket2 & BasisType::bitmask(i));
						if (s1i>0) s1i=1;
						if (s2i>0) s2i=1;
						
						// Hubbard term
						if (s1i>0 && s2i>0 ) s += mp_.hubbardU[i];
								
						// Potential term
						if (s1i>0) s += mp_.potentialV[i];
						if (s2i>0) s += mp_.potentialV[i];
						
						// Hopping term (only count how many non-zero)
						for (size_t j=0;j<nsite;j++) {
							if (j<i) continue;
							RealType tmp = hoppings(i,j);
							if (tmp==0) continue;
							
							WordType s1j= (ket1 & BasisType::bitmask(j));
							WordType s2j= (ket2 & BasisType::bitmask(j));
							if (s1j>0) s1j=1;
							if (s2j>0) s2j=1;
							if (s1i+s1j==1) nzero++;
							if (s2i+s2j==1) nzero++;
						}
					}
					// cout<<"diag of ("<<ispace1<<","<<ispace2<<"): "<<s<<endl;
					diag(ispace2,ispace1)=s;
					nzero++;
				}
			}

			nzero++;
			return nzero;
		}
		

		size_t perfectIndex(const BasisType& basis1,const BasisType& basis2,WordType ket1,WordType ket2) const
		{
			size_t hilbert2=basis2.size();
			size_t n1 = basis1.perfectIndex(ket1);
			size_t n2 = basis2.perfectIndex(ket2);

			return n2 + n1*hilbert2;
		}
		

		int doSign(WordType a, WordType b,size_t i,size_t j,size_t sector) const
		{
			if (i > j) {
				std::cerr<<"FATAL: At doSign\n";
				std::cerr<<"INFO: i="<<i<<" j="<<j<<std::endl;
				std::cerr<<"AT: "<<__FILE__<<" : "<<__LINE__<<std::endl;
				throw std::runtime_error("HubbardOneOrbital::doSign(...)\n");
			}

			if (sector==SPIN_UP) {
				WordType mask = a;
				mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
				int s=(PsimagLite::BitManip::count(mask) & 1) ? -1 : 1; // Parity of up between i and j 
				// Is there an up at i?(killed due to Hermitian)
				// if (BasisType::bitmask(i) & a) s = -s;
				return s;
			}
			WordType mask = b;
			mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
			int s=(PsimagLite::BitManip::count(mask) & 1) ? -1 : 1; // Parity of down between i and j 
			//if (sector==SPIN_DOWN) {
			
			// Is there a down at i? (killed due to Hermitian)
			// if (BasisType::bitmask(i) & b) s = -s;
			//}

			return s;
		}
		

		int doSign(WordType a, WordType b,size_t ind,size_t sector) const
		{
			if (sector==SPIN_UP) {
				if (ind==0) return 1;
	
				// ind>0 from now on
				size_t i = 0;
				size_t j = ind;
				WordType mask = a;
				mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
				int s=(PsimagLite::BitManip::count(mask) & 1) ? -1 : 1; // Parity of up between i and j
				return s;
			}
			int s=(PsimagLite::BitManip::count(a) & 1) ? -1 : 1; // Parity of up
			if (ind==0) return s;
	
			// ind>0 from now on
			size_t i = 0;
			size_t j = ind;
			WordType mask = b;
			mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
			s=(PsimagLite::BitManip::count(mask) & 1) ? -1 : 1; // Parity of up between i and j
			return s;
			//if (sector==SPIN_DOWN) {
			// is there an up at j?
			//if (BasisType::bitmask(j) & a) s = -s;
			//}

			return s;
		}
		
		//! Gf. related functions below:

		std::pair<int,int> getNewParts(size_t type,size_t spin) const
		{
			throw std::runtime_error("getNewParts: unimplemented");
		}
		
		
		void getModifiedStates(
				std::vector<RealType>& modifVector,
				size_t type,
				size_t isite,
				size_t jsite,
				size_t spin) const
		{
			throw std::runtime_error("getModifiedStates: unimplemented");
		}


		template<typename ContinuedFractionType>
		void calcGf(
				ContinuedFractionType& cf,
				const std::pair<int,int>& newParts,
				const std::vector<RealType>& modifVector,
				size_t type,
				size_t spin) const
		{
			throw std::runtime_error("calcGf: unimplemented");
		}

		const ParametersType& mp_;
		const GeometryType& geometry_;
		BasisType basis1_;
		BasisType basis2_;
		
	}; // class HubbardOneOrbital 
} // namespace LanczosPlusPlus
#endif

