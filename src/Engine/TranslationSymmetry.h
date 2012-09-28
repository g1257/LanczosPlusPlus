
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
#include "SparseVector.h"

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

	const RealType& operator()(size_t i,size_t j) const
	{
		return data_(i,j);
	}

private:

	PsimagLite::Matrix<RealType> data_;
};

template<typename GeometryType,typename BasisType,typename KspaceType>
class ClassRepresentatives {

	typedef typename BasisType::WordType WordType;
	typedef TranslationItem ItemType;

public:

	ClassRepresentatives(const BasisType& basis,const GeometryType& geometry,const KspaceType& kspace)
		: basis_(basis),geometry_(geometry),data_(basis.size())
	{
		size_t hilbert = basis.size();
		std::vector<bool> seen(hilbert,false);

		for (size_t ispace=0;ispace<hilbert;ispace++) {

			for (size_t k=0;k<kspace.size();k++) {
				std::vector<WordType> y = translateInternal(ispace,k);
				size_t yIndex = basis.perfectIndex(y);
				if (!seen[yIndex]) {
					seen[yIndex]=true;
					data_[yIndex].type = k;
					data_[yIndex].i = ispace;
				}
			}
		}

		for (size_t i=0;i<data_.size();i++) {
			if (data_[i].type!=0) continue;
			reps_.push_back(i);
			std::cout<<"Rep "<<basis(i,0)<<" "<<basis(i,1)<<"\n";
		}
		std::cout<<"Total Representatives="<<reps_.size()<<"\n";
	}

	size_t size() const { return reps_.size(); }

	const size_t& operator[](size_t i) const { return reps_[i]; }

	size_t translate(size_t state,size_t k) const
	{
		for (size_t i=0;i<data_.size();i++) {
			if (data_[i].i == state && data_[i].type == k) return i;
		}
//		assert(false);
//		throw std::runtime_error("TranslationSymmetry: translate\n");
		return state;
	}

private:

	std::vector<WordType> translateInternal(size_t state,size_t k) const
	{
		size_t numberOfDofs = basis_.dofs();
		std::vector<WordType> y(numberOfDofs);

		for (size_t dof=0;dof<numberOfDofs;dof++) {
			WordType x = basis_(state,dof);
			y[dof] = translateInternal2(x,k);
		}
		std::cout<<"translation of "<<basis_(state,0)<<" "<<basis_(state,1);
		std::cout<<"is "<<y[0]<<" "<<y[1]<<"\n";
		return y;
	}

	WordType translateInternal2(WordType state,size_t k) const
	{
		size_t numberOfSites = geometry_.numberOfSites();
		size_t termId = 0;
		WordType x = state;
		WordType y = 0;
		size_t diry = 1;
		for (size_t site=0;site<numberOfSites;site++) {
			size_t tSite = geometry_.translate(site,diry,k,termId);
			size_t thisSiteContent = x & 1;
			x >>=1; // go to next site
			addTo(y,thisSiteContent,tSite);
			if (!x) break;
		}

		return y;
	}

	void addTo(WordType& yy,size_t what,size_t site) const
	{
		if (what==0) return;
		WordType mask = (1<<site);
		yy |= mask;
	}

	const BasisType& basis_;
	const GeometryType& geometry_;
	std::vector<TranslationItem> data_;
	std::vector<size_t> reps_;
};

	template<typename GeometryType,typename BasisType>
	class TranslationSymmetry  {

		typedef typename GeometryType::RealType RealType;
		typedef std::complex<RealType> ComplexType;
		typedef typename BasisType::WordType WordType;
		typedef PsimagLite::SparseVector<ComplexType> SparseVectorType;
		typedef Kspace<RealType> KspaceType;
		typedef ClassRepresentatives<GeometryType,BasisType,KspaceType> ClassRepresentativesType;
		typedef std::pair<std::vector<size_t> ,size_t> BufferItemType;

	public:

		typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
		typedef std::vector<ComplexType> VectorType;

		TranslationSymmetry(const BasisType& basis,const GeometryType& geometry)
		: progress_("TranslationSymmetry",0),
		  transform_(basis.size(),basis.size()),
		  kspace_(geometry.length(1,0)),
		  matrixStored_(kspace_.size()),
		  pointer_(0)
		{
			ClassRepresentativesType reps(basis,geometry,kspace_);

			size_t hilbert = basis.size();
			std::vector<SparseVectorType> bag;
			for (size_t k=0;k<kspace_.size();k++) {
				for (size_t ispace=0;ispace<hilbert;ispace++) {
					std::vector<ComplexType> v(hilbert);
					eikrTr(v,ispace,k,reps);
					SparseVectorType sparseV(v);
					sparseV.sort();
					if (!checkForOrthogonality(sparseV,bag)) continue;
					bag.push_back(sparseV);
				}
			}

			setTransform(bag);
//			checkTransform();
		}

		template<typename SomeModelType>
		void init(const SomeModelType& model,const BasisType& basis)
		{
			PsimagLite::CrsMatrix<RealType> matrix2;
			model.setupHamiltonian(matrix2,basis);
			transformMatrix(matrixStored_,matrix2);
		}

		size_t rank() const { return matrixStored_[pointer_].row(); }

		template<typename SomeVectorType>
		void matrixVectorProduct(SomeVectorType &x, SomeVectorType const &y) const
		{
			return matrixStored_[pointer_].matrixVectorProduct(x,y);
		}

		void transformMatrix(std::vector<SparseMatrixType>& matrix1,const PsimagLite::CrsMatrix<RealType>& matrix) const
		{
			SparseMatrixType rT;
			transposeConjugate(rT,transform_);
			
			printFullMatrix(matrix,"originalHam");
			SparseMatrixType tmp;
			multiply(tmp,matrix,rT);

			SparseMatrixType matrix2;
			multiply(matrix2,transform_,tmp);

			if (matrix2.row()<40)
				printFullMatrix(matrix2,"HamiltonianTransformed");
			matrix1.clear();
			split(matrix1,matrix2);
			assert(matrix1.size()==kspace_.size());
		}

		void transformGs(VectorType& gs,size_t offset)
		{
			VectorType gstmp(transform_.row(),0);

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

		size_t sectors() const { return kspace_.size(); }

		void setPointer(size_t p) { pointer_=p; }

		std::string name() const { return "translation"; }

	private:

		void addTo(WordType& yy,size_t what,size_t site) const
		{
			if (what==0) return;
			WordType mask = (1<<site);
			yy |= mask;
		}

		void setTransform(const std::vector<SparseVectorType>& bag)
		{
			size_t counter = 0;
			for (size_t row=0;row<bag.size();row++) {
				transform_.setRow(row,counter);
				for (size_t k=0;k<bag[row].indices();k++) {
					transform_.pushCol(bag[row].index(k));
					transform_.pushValue(bag[row].value(k));
					counter++;
				}
			}

			transform_.setRow(transform_.row(),counter);
			transform_.checkValidity();
			if (transform_.row()<40)
				printFullMatrix(transform_,"transform");
		}

		void eikrTr(std::vector<ComplexType>& v,size_t ispace,size_t k,const ClassRepresentativesType& reps) const
		{
			for (size_t r=0;r<kspace_.size();r++) {
				size_t jspace = reps.translate(ispace,k);
				RealType tmp = 2*M_PI*k*r/RealType(kspace_.size());
				v[jspace] = ComplexType(cos(tmp),sin(tmp));
			}
		}

		bool checkForOrthogonality(const SparseVectorType& sparseV,const std::vector<SparseVectorType>& bag) const
		{
			for (size_t i=0;i<bag.size();i++) {
				ComplexType sp = sparseV.scalarProduct(bag[i]);

				if (std::norm(sp)>1e-8) return false;
			}
			return true;
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

		void split(std::vector<SparseMatrixType>& matrix,const SparseMatrixType& matrix2) const
		{
			size_t blockSize = matrix2.row()/kspace_.size();
			std::cout<<"BLOCKSIZE="<<blockSize<<"\n";
			for (size_t i=0;i<kspace_.size();i++) {
				SparseMatrixType m(blockSize,blockSize);
				size_t offset = i*blockSize;
				size_t counter = 0;
				for (size_t row=0;row<blockSize;row++) {
					m.setRow(row,counter);
					size_t globalRow = row + offset;
					for (int k=matrix2.getRowPtr(globalRow);k<matrix2.getRowPtr(globalRow+1);k++) {
						ComplexType val = matrix2.getValue(k);
						if (std::norm(val)<1e-8) continue;
						size_t globalCol = matrix2.getCol(k);
//						if (globalCol<offset) continue; // <<-- FIXME
						assert(globalCol>=offset);
						size_t col = globalCol - offset;
//						if (col>=blockSize) continue; // <<-- FIXME
						assert(col<blockSize);
						m.pushCol(col);
						m.pushValue(val);
						counter++;
					}
				}
				m.setRow(blockSize,counter);
				m.checkValidity();
				matrix.push_back(m);
			}
		}
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

		PsimagLite::ProgressIndicator progress_;
		SparseMatrixType transform_;
		KspaceType kspace_;
		std::vector<SparseMatrixType> matrixStored_;
		size_t pointer_;
//		SparseMatrixType s_;
	}; // class TranslationSymmetry
} // namespace Dmrg

#endif  // TRANSLATION_SYMM_H
