
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

#ifndef REFLECTION_SYMM_H
#define REFLECTION_SYMM_H
#include <iostream>
#include "ProgressIndicator.h"
#include "CrsMatrix.h"
#include "Vector.h"

namespace LanczosPlusPlus {

	class ReflectionItem {

	public:

		enum { DIAGONAL,PLUS,MINUS};

		ReflectionItem(size_t ii)
		: i(ii),j(ii),type(DIAGONAL)
		{}

		ReflectionItem(size_t ii,size_t jj,size_t type1)
		: i(ii),j(jj),type(type1)
		{}

		size_t i,j,type;

	}; // class ReflectionItem

	bool operator==(const ReflectionItem& item1,const ReflectionItem& item2)
	{
		if (item1.type!=item2.type) return false;

		if (item1.i==item2.j && item1.j==item2.i) return true;

		return (item1.i==item2.i && item1.j==item2.j);
	}

	template<typename GeometryType,typename BasisType>
	class ReflectionSymmetry  {

		typedef typename GeometryType::RealType RealType;
		typedef typename BasisType::WordType WordType;
		typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;


		typedef ReflectionItem ItemType;

	public:

		ReflectionSymmetry(const BasisType& basis,const GeometryType& geometry)
		: progress_("ReflectionSymmetry",0),
		  transform_(basis.size(),basis.size())
//		  s_(basis.size(),basis.size()) // needed only for debugging
		{
			size_t hilbert = basis.size();
			size_t numberOfDofs = basis.dofs();
			size_t numberOfSites = geometry.numberOfSites();
			size_t termId = 0;
//			size_t counter=0;
			std::vector<ItemType> buffer;
			for (size_t ispace=0;ispace<hilbert;ispace++) {
				std::vector<WordType> y(numberOfDofs,0);
				for (size_t dof=0;dof<numberOfDofs;dof++) {
					WordType x = basis(ispace,dof);
					for (size_t site=0;site<numberOfSites;site++) {
						size_t reflectedSite = geometry.findReflection(site,termId);
						size_t thisSiteContent = x & 1;
						x >>=1; // go to next site
						addTo(y[dof],thisSiteContent,reflectedSite);
						if (!x) break;
					}
				}

				size_t yIndex = basis.perfectIndex(y);
//				s_.setRow(ispace,counter);
//				s_.pushCol(yIndex);
//				s_.pushValue(1.0);
//				counter++;
				if (yIndex==ispace) { // then S|psi> = |psi>
					ItemType item1(ispace);
					buffer.push_back(item1);
					continue;
				}
				// S|psi> != |psi>
				// Add normalized +
				ItemType item2(ispace,yIndex,ItemType::PLUS);
				buffer.push_back(item2);

				// Add normalized -
				ItemType item3(ispace,yIndex,ItemType::MINUS);
				buffer.push_back(item3);
			}
//			s_.setRow(s_.rank(),counter);
			setTransform(buffer);
//			checkTransform();
		}

		void transform(SparseMatrixType& matrix2,const SparseMatrixType& matrix) const
		{
//			PsimagLite::Matrix<RealType> fullMatrix;
//			crsMatrixToFullMatrix(fullMatrix,matrix);
//			std::cerr<<"-----------\n";
//			std::cerr<<fullMatrix<<"\n";

//			ReflectionSymmetryType rs(basis_,geometry_);

//			const SparseMatrixType& r = rs();
			SparseMatrixType rT;
			transposeConjugate(rT,transform_);

			SparseMatrixType tmp;
			multiply(tmp,matrix,rT);
//			PsimagLite::Matrix<RealType> mtmp;
//			crsMatrixToFullMatrix(mtmp,tmp);
//			std::cerr<<"-----------\n";
//			std::cerr<<mtmp;

			multiply(matrix2,transform_,tmp);
//			PsimagLite::Matrix<RealType> mtmp2;
//			crsMatrixToFullMatrix(mtmp2,tmp2);
//			std::cerr<<"-----------\n";
//			std::cerr<<mtmp2;

//			throw std::runtime_error("testing\n");
		}

		//const SparseMatrixType operator()() const { return transform_; }

//		const SparseMatrixType reflectionSymmetry() const { return s_; }

	private:

		void checkTransform() const
		{
			SparseMatrixType transformTc;
			transposeConjugate(transformTc,transform_);

//			SparseMatrixType tmp;
//			multiply(tmp,s_,s_);
			PsimagLite::Matrix<RealType> mtmp;
			crsMatrixToFullMatrix(mtmp,transform_);
			std::cerr<<"&&&&&&&&&&&&&&&&&&&&&&&\n";
			std::cerr<<mtmp;
//			throw std::runtime_error("checking\n");
		}

		void addTo(WordType& yy,size_t what,size_t site) const
		{
			if (what==0) return;
			WordType mask = (1<<site);
			yy |= mask;
		}

		void setTransform(const std::vector<ItemType>& buffer2)
		{
			std::vector<ItemType> buffer;
			makeUnique(buffer,buffer2);
			assert(buffer.size()==transform_.rank());
			size_t counter = 0;
			RealType oneOverSqrt2 = 1.0/sqrt(2.0);
			RealType sign = 1.0;
			size_t row = 0;
			for (size_t i=0;i<buffer.size();i++) {
				if (buffer[i].type==ItemType::MINUS) continue;
				transform_.setRow(row++,counter);
				switch(buffer[i].type) {
				case ItemType::DIAGONAL:
					transform_.pushCol(buffer[i].i);
					transform_.pushValue(1);
					counter++;
					break;
				case ItemType::PLUS:
					transform_.pushCol(buffer[i].i);
					transform_.pushValue(oneOverSqrt2);
					counter++;
					transform_.pushCol(buffer[i].j);
					transform_.pushValue(oneOverSqrt2);
					counter++;
					break;
				}
			}

			for (size_t i=0;i<buffer.size();i++) {
				if (buffer[i].type!=ItemType::MINUS) continue;
				transform_.setRow(row++,counter);
				transform_.pushCol(buffer[i].i);
				transform_.pushValue(oneOverSqrt2*sign);
				counter++;
				transform_.pushCol(buffer[i].j);
				transform_.pushValue(-oneOverSqrt2*sign);
				counter++;
			}
			transform_.setRow(transform_.rank(),counter);
		}

		void makeUnique(std::vector<ItemType>& dest,const std::vector<ItemType>& src)
		{
			size_t zeros=0;
			size_t pluses=0;
			size_t minuses=0;
			for (size_t i=0;i<src.size();i++) {
				ItemType item = src[i];
				int x =  PsimagLite::isInVector(dest,item);
				if (x>=0) continue;
//				if (item.type ==ItemType::PLUS) {
//					size_t i = item.i;
//					size_t j = item.j;
//					ItemType item2(j,i,ItemType::PLUS);
//					x = PsimagLite::isInVector(dest,item2);
//					if (x>=0) continue;
//				}
				if (item.type==ItemType::DIAGONAL) zeros++;
				if (item.type==ItemType::PLUS) pluses++;
				if (item.type==ItemType::MINUS) minuses++;

				dest.push_back(item);
			}
			std::ostringstream msg;
			msg<<pluses<<" +, "<<minuses<<" -, "<<zeros<<" zeros.";
			progress_.printline(msg,std::cout);
		}

		void isIdentity(const SparseMatrixType& s,const std::string& label) const
		{
			std::cerr<<"Checking label="<<label<<"\n";
			for (size_t i=0;i<s.rank();i++) {
				for (int k=s.getRowPtr(i);k<s.getRowPtr(i+1);k++) {
					size_t col = s.getCol(k);
					RealType val = s.getValue(k);
					if (col==i) assert(isAlmostZero(val-1.0));
					else assert(isAlmostZero(val));
				}
			}
		}

		bool isAlmostZero(const RealType& x) const
		{
			return (fabs(x)<1e-6);
		}

		PsimagLite::ProgressIndicator progress_;
		SparseMatrixType transform_;
//		SparseMatrixType s_;
	}; // class ReflectionSymmetry
} // namespace Dmrg

#endif  // REFLECTION_SYMM_H
