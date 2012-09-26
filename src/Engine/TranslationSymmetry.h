
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

#ifndef TRANSLATION_SYMM_H
#define TRANSLATION_SYMM_H
#include <iostream>
#include "ProgressIndicator.h"
#include "CrsMatrix.h"
#include "Vector.h"

namespace LanczosPlusPlus {

struct TranslationItem {
	size_t type;
	size_t i;
};

template<typename RealType>
class Kspace {

public:

	Kspace(size_t len) : data_(len,len)
	{

		for (size_t i=0;i<len-1;i++) {
			data_(i,i+1) = 1;
		}
		data_(len-1,0) =1;
		std::vector<RealType> e(len);
		diag(data_,e,'V');
	}

	size_t size() const { return data_.n_row(); }

private:

	PsimagLite::Matrix<RealType> data_;
};

template<typename GeometryType,typename BasisType,typename KspaceType>
class ClassRepresentatives {

	typedef typename BasisType::WordType WordType;
	typedef TranslationItem ItemType;

public:

	ClassRepresentatives(const BasisType& basis,const GeometryType& geometry,const KspaceType& kspace)
		: data_(basis.size())
	{
//		size_t hilbert = basis.size();
//		std::vector<bool> seen(hilbert,false);

//		for (size_t ispace=0;ispace<hilbert;ispace++) {

//			for (size_t k=0;k<kspace.size();k++) {
//				std::vector<WordType> y = kspace.translate(ispace,k);
//				size_t yIndex = basis.perfectIndex(y);
//				if (!seen[yIndex]) {
//					seen[yIndex]=true;
//					data_[yIndex].type = k;
//					data_[yIndex].i = ispace;
//				}
//			}
//		}

//		for (size_t i=0;i<data_.size();i++) {
//			if (data_[i].type!=0) continue;
//			reps_.push_back(i);
//		}
		throw std::runtime_error("TranslationSymmetry: more work needed\n");
	}

	size_t size() const { return reps_.size(); }

	const size_t& operator[](size_t i) const { return reps_[i]; }

	size_t translate(size_t state,size_t k) const
	{
		for (size_t i=0;i<data_.size();i++) {
			if (data_[i].i == state && data_[i].type == k) return i;
		}
		assert(false);
		throw std::runtime_error("TranslationSymmetry: translate\n");
	}

private:

	std::vector<TranslationItem> data_;
	std::vector<size_t> reps_;
};

	template<typename GeometryType,typename BasisType>
	class TranslationSymmetry  {

		typedef typename GeometryType::RealType RealType;
		typedef typename BasisType::WordType WordType;
		typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;
		typedef std::vector<RealType> VectorType;

		typedef Kspace<RealType> KspaceType;
		typedef ClassRepresentatives<GeometryType,BasisType,KspaceType> ClassRepresentativesType;
		typedef std::pair<std::vector<size_t> ,size_t> BufferItemType;

	public:

		TranslationSymmetry(bool enabled,const BasisType& basis,const GeometryType& geometry)
		: enabled_(enabled),
		  progress_("TranslationSymmetry",0),
		  transform_(basis.size(),basis.size()),
		  plusSector_(0)
		{
			if (!enabled_) return;

			size_t termId = 0;
			KspaceType kspace(geometry.length(1,termId));
			ClassRepresentativesType classReps(basis,geometry,kspace);

			std::vector<BufferItemType> buffer;
			for (size_t st=0;st<classReps.size();st++) {

				std::vector<size_t> translatedIndices(kspace.size());
				for (size_t i=0;i<kspace.size();i++) {
					translatedIndices[i]  = classReps.translate(classReps[st],i);
				}
				for (size_t i=0;i<kspace.size();i++) {
					BufferItemType item;
					item.first = translatedIndices;
					item.second = i;
					buffer.push_back(item);
				}
			}

			setTransform(buffer,kspace);
//			checkTransform();
		}

		void transformMatrix(std::vector<SparseMatrixType>& matrix1,const SparseMatrixType& matrix) const
		{
			if (!enabled_) {
				throw std::runtime_error("TranslationSymmetry: transform(...) called on disabled\n");
			}
			SparseMatrixType rT;
			transposeConjugate(rT,transform_);
			
			printFullMatrix(matrix,"originalHam");
			SparseMatrixType tmp;
			multiply(tmp,matrix,rT);

			SparseMatrixType matrix2;
			multiply(matrix2,transform_,tmp);

			assert(matrix1.size()==2);
			//split(matrix1[0],matrix1[1],matrix2);
			throw std::runtime_error("TranslationSymmetry: needs more work\n");
		}

		void transformGs(VectorType& gs,size_t offset)
		{
			if (!enabled_) return;

			std::vector<RealType> gstmp(transform_.row(),0);

			for (size_t i=0;i<gs.size();i++) {
				assert(i+offset<gstmp.size());
				gstmp[i+offset]=gs[i];
			}
			SparseMatrixType rT;
			transposeConjugate(rT,transform_);
			gs.clear();
			gs.resize(transform_.row());
			multiply(gs,rT,gstmp);
		}

		size_t sectors() const { return 2; }

	private:

		void addTo(WordType& yy,size_t what,size_t site) const
		{
			if (what==0) return;
			WordType mask = (1<<site);
			yy |= mask;
		}

		void setTransform(const std::vector<BufferItemType>& buffer,const KspaceType& kspace)
		{
			throw std::runtime_error("TranslationSymmetry: needs more work\n");

//			assert(buffer.size()==transform_.row());
//			size_t counter = 0;

//			size_t row = 0;

//			for (size_t k=0;k<kspace.size();k++) {
//				for (size_t i=0;i<buffer.size();i++) {
//					if (buffer[i].second==k) continue;

//					transform_.setRow(row++,counter);

//					transform_.pushCol(buffer[i].first);
//					transform_.pushValue(kspace(k,i));
//					counter++;
//				}
//			}

//			transform_.setRow(transform_.row(),counter);
//			transform_.checkValidity();
		}

//		void makeUnique(std::vector<ItemType>& dest,const std::vector<ItemType>& src)
//		{
//			size_t zeros=0;
//			size_t pluses=0;
//			size_t minuses=0;
//			for (size_t i=0;i<src.size();i++) {
//				ItemType item = src[i];
//				int x =  PsimagLite::isInVector(dest,item);
//				if (x>=0) continue;

//				if (item.type==ItemType::DIAGONAL) zeros++;
//				if (item.type==ItemType::PLUS) pluses++;
//				if (item.type==ItemType::MINUS) minuses++;

//				dest.push_back(item);
//			}
//			std::ostringstream msg;
//			msg<<pluses<<" +, "<<minuses<<" -, "<<zeros<<" zeros.";
//			progress_.printline(msg,std::cout);
//			plusSector_ = zeros + pluses;
//		}

//		void isIdentity(const SparseMatrixType& s,const std::string& label) const
//		{
//			std::cerr<<"Checking label="<<label<<"\n";
//			for (size_t i=0;i<s.rank();i++) {
//				for (int k=s.getRowPtr(i);k<s.getRowPtr(i+1);k++) {
//					size_t col = s.getCol(k);
//					RealType val = s.getValue(k);
//					if (col==i) assert(isAlmostZero(val-1.0));
//					else assert(isAlmostZero(val));
//				}
//			}
//		}

//		bool isAlmostZero(const RealType& x) const
//		{
//			return (fabs(x)<1e-6);
//		}

//		void split(SparseMatrixType& matrixA,SparseMatrixType& matrixB,const SparseMatrixType& matrix) const
//		{
//			size_t counter = 0;
//			matrixA.resize(plusSector_,plusSector_);
//			for (size_t i=0;i<plusSector_;i++) {
//				matrixA.setRow(i,counter);
//				for (int k=matrix.getRowPtr(i);k<matrix.getRowPtr(i+1);k++) {
//					size_t col = matrix.getCol(k);
//					RealType val = matrix.getValue(k);
//					if (col<plusSector_) {
//						matrixA.pushCol(col);
//						matrixA.pushValue(val);
//						counter++;
//						continue;
//					}
//					if (!isAlmostZero(val)) {
//						std::string s(__FILE__);
//						s += " Hamiltonian has no reflection symmetry.";
//						throw std::runtime_error(s.c_str());
//					}
//				}
//			}
//			matrixA.setRow(plusSector_,counter);

//			size_t rank = matrix.row();
//			size_t minusSector=rank-plusSector_;
//			matrixB.resize(minusSector,minusSector);
//			counter=0;
//			for (size_t i=plusSector_;i<rank;i++) {
//				matrixB.setRow(i-plusSector_,counter);
//				for (int k=matrix.getRowPtr(i);k<matrix.getRowPtr(i+1);k++) {
//					size_t col = matrix.getCol(k);
//					RealType val = matrix.getValue(k);
//					if (col>=plusSector_) {
//						matrixB.pushCol(col-plusSector_);
//						matrixB.pushValue(val);
//						counter++;
//						continue;
//					}
//					if (!isAlmostZero(val)) {
//						std::string s(__FILE__);
//						s += " Hamiltonian has no reflection symmetry.";
//						throw std::runtime_error(s.c_str());
//					}
//				}
//			}
//			matrixB.setRow(minusSector,counter);
//		}

		bool enabled_;
		PsimagLite::ProgressIndicator progress_;
		SparseMatrixType transform_;
		size_t plusSector_;
//		SparseMatrixType s_;
	}; // class TranslationSymmetry
} // namespace Dmrg

#endif  // TRANSLATION_SYMM_H