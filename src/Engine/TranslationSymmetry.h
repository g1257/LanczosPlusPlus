
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
	SizeType type;
	SizeType i;
};

template<typename RealType>
class Kspace {

public:

	Kspace(SizeType len) : data_(len,len),blockSizes_(len,0)
	{

		for (SizeType i=0;i<len-1;i++) {
			data_(i,i+1) = 1;
		}
		data_(len-1,0) =1;
		typename PsimagLite::Vector<RealType>::Type e(len);
		diag(data_,e,'V');
	}

	SizeType size() const { return data_.n_row(); }

	const RealType& operator()(SizeType i,SizeType j) const
	{
		return data_(i,j);
	}

	void setBlockSize(SizeType k,SizeType s)
	{
		std::cout<<"BLOCKSZIZE["<<k<<"]="<<s<<"\n";
		blockSizes_[k]=s;
	}

	SizeType blockSizes(SizeType k) const
	{
		return blockSizes_[k];
	}

	SizeType blockSize() const
	{
		SizeType sum = 0;
		for (SizeType i=0;i<blockSizes_.size();i++)
			sum += blockSizes_[i];
		return sum;
	}

private:

	PsimagLite::Matrix<RealType> data_;
	PsimagLite::Vector<SizeType>::Type blockSizes_;
};

template<typename GeometryType,typename BasisType,typename KspaceType>
class ClassRepresentatives {

	typedef typename BasisType::WordType WordType;
	typedef TranslationItem ItemType;

public:

	ClassRepresentatives(const BasisType& basis,const GeometryType& geometry,const KspaceType& kspace)
		: basis_(basis),geometry_(geometry),data_(basis.size())
	{
		SizeType hilbert = basis.size();
		PsimagLite::Vector<bool>::Type seen(hilbert,false);

		for (SizeType ispace=0;ispace<hilbert;ispace++) {

			for (SizeType k=0;k<kspace.size();k++) {
				typename PsimagLite::Vector<WordType>::Type y = translateInternal(ispace,k);
				SizeType yIndex = basis.perfectIndex(y);
				if (!seen[yIndex]) {
					seen[yIndex]=true;
					data_[yIndex].type = k;
					data_[yIndex].i = ispace;
				}
			}
		}

//		for (SizeType i=0;i<data_.size();i++) {
//			if (data_[i].type!=0) continue;
//			reps_.push_back(i);
//			std::cout<<"Rep "<<basis(i,0)<<" "<<basis(i,1)<<"\n";
//		}
//		std::cout<<"Total Representatives="<<reps_.size()<<"\n";
	}

//	SizeType size() const { return reps_.size(); }

//	const SizeType& operator[](SizeType i) const { return reps_[i]; }

	SizeType translate(SizeType state,SizeType k) const
	{
		for (SizeType i=0;i<data_.size();i++) {
			if (data_[i].i == state && data_[i].type == k) return i;
		}
//		assert(false);
//		throw std::runtime_error("TranslationSymmetry: translate\n");
		return state;
	}

private:

	typename PsimagLite::Vector<WordType>::Type translateInternal(SizeType state,SizeType k) const
	{
		SizeType numberOfDofs = basis_.dofs();
		typename PsimagLite::Vector<WordType>::Type y(numberOfDofs);

		for (SizeType dof=0;dof<numberOfDofs;dof++) {
			WordType x = basis_(state,dof);
			y[dof] = translateInternal2(x,k);
		}
		std::cout<<"translation of "<<basis_(state,0)<<" "<<basis_(state,1);
		std::cout<<"is "<<y[0]<<" "<<y[1]<<"\n";
		return y;
	}

	WordType translateInternal2(WordType state,SizeType k) const
	{
		SizeType numberOfSites = geometry_.numberOfSites();
		SizeType termId = 0;
		WordType x = state;
		WordType y = 0;
		SizeType diry = 1;
		for (SizeType site=0;site<numberOfSites;site++) {
			SizeType tSite = geometry_.translate(site,diry,k,termId);
			SizeType thisSiteContent = x & 1;
			x >>=1; // go to next site
			addTo(y,thisSiteContent,tSite);
			if (!x) break;
		}

		return y;
	}

	void addTo(WordType& yy,SizeType what,SizeType site) const
	{
		if (what==0) return;
		WordType mask = (1<<site);
		yy |= mask;
	}

	const BasisType& basis_;
	const GeometryType& geometry_;
	PsimagLite::Vector<TranslationItem>::Type data_;
//	PsimagLite::Vector<SizeType>::Type reps_;
};

	template<typename GeometryType,typename BasisType>
	class TranslationSymmetry  {

		typedef typename GeometryType::ComplexOrRealType ComplexOrRealType;
		typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
		typedef typename BasisType::WordType WordType;
		typedef PsimagLite::SparseVector<ComplexOrRealType> SparseVectorType;
		typedef Kspace<RealType> KspaceType;
		typedef ClassRepresentatives<GeometryType,BasisType,KspaceType> ClassRepresentativesType;
		typedef std::pair<PsimagLite::Vector<SizeType>::Type ,SizeType> BufferItemType;

	public:

		typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
		typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

		TranslationSymmetry(const BasisType& basis,const GeometryType& geometry)
		: progress_("TranslationSymmetry"),
		  transform_(basis.size(),basis.size()),
		  kspace_(geometry.length(1,0)),
		  matrixStored_(kspace_.size()),
		  pointer_(0)
		{
			ClassRepresentativesType reps(basis,geometry,kspace_);

			SizeType hilbert = basis.size();
			typename PsimagLite::Vector<SparseVectorType>::Type bag;
			for (SizeType k=0;k<kspace_.size();k++) {
				SizeType blockSize = 0;
				for (SizeType ispace=0;ispace<hilbert;ispace++) {
					typename PsimagLite::Vector<ComplexOrRealType>::Type v(hilbert);
					eikrTr(v,ispace,k,reps);
					SparseVectorType sparseV(v);
					sparseV.sort();
					if (!checkForOrthogonality(sparseV,bag)) continue;
					bag.push_back(sparseV);
					blockSize++;
				}
				kspace_.setBlockSize(k,blockSize);
			}

			if (kspace_.blockSize()!=hilbert) {
				std::cout<<"Blocksizes summed="<<kspace_.blockSize()<<" but hilbert="<<hilbert<<"\n";
				throw std::runtime_error("error!\n");
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

		SizeType rank() const { return matrixStored_[pointer_].row(); }

		template<typename SomeVectorType>
		void matrixVectorProduct(SomeVectorType &x, SomeVectorType const &y) const
		{
			return matrixStored_[pointer_].matrixVectorProduct(x,y);
		}

		void transformMatrix(typename PsimagLite::Vector<SparseMatrixType>::Type& matrix1,
		                     const PsimagLite::CrsMatrix<RealType>& matrix) const
		{
			SparseMatrixType rT;
			transposeConjugate(rT,transform_);
			
			if (matrix.row()<40) printFullMatrix(matrix,"originalHam");
			SparseMatrixType tmp;
			multiply(tmp,matrix,rT);

			SparseMatrixType matrix2;
			multiply(matrix2,transform_,tmp);

			if (matrix2.row()<40)
				printFullMatrix(matrix2,"HamiltonianTransformed");
			matrix1.clear();
			split(matrix1,matrix2);
//			assert(matrix1.size()==kspace_.size());
		}

		void transformGs(VectorType& gs,SizeType offset)
		{
			VectorType gstmp(transform_.row(),0);

			for (SizeType i=0;i<gs.size();i++) {
				assert(i+offset<gstmp.size());
				gstmp[i+offset]=gs[i];
			}
			SparseMatrixType rT;
			transposeConjugate(rT,transform_);
			gs.clear();
			gs.resize(transform_.row());
			multiply(gs,rT,gstmp);
		}

		SizeType sectors() const { return kspace_.size(); }

		void setPointer(SizeType p) { pointer_=p; }

		PsimagLite::String name() const { return "translation"; }

	private:

		void addTo(WordType& yy,SizeType what,SizeType site) const
		{
			if (what==0) return;
			WordType mask = (1<<site);
			yy |= mask;
		}

		void setTransform(const typename PsimagLite::Vector<SparseVectorType>::Type& bag)
		{
			SizeType counter = 0;
			for (SizeType row=0;row<bag.size();row++) {
				transform_.setRow(row,counter);
				for (SizeType k=0;k<bag[row].indices();k++) {
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

		void eikrTr(typename PsimagLite::Vector<std::complex<RealType> >::Type& v,
		            SizeType ispace,
		            SizeType k,
		            const ClassRepresentativesType& reps) const
		{
			for (SizeType r=0;r<kspace_.size();r++) {
				SizeType jspace = reps.translate(ispace,r);
				RealType tmp = 2*M_PI*k*r/RealType(kspace_.size());
				v[jspace] = ComplexOrRealType(cos(tmp),sin(tmp));
			}
		}

		void eikrTr(typename PsimagLite::Vector<RealType>::Type& v,
		            SizeType ispace,
		            SizeType k,
		            const ClassRepresentativesType& reps) const
		{
			throw PsimagLite::RuntimeError("eikrTr: not for real template\n");
		}

		bool checkForOrthogonality(const SparseVectorType& sparseV,
		                           const typename PsimagLite::Vector<SparseVectorType>::Type& bag) const
		{
			for (SizeType i=0;i<bag.size();i++) {
				ComplexOrRealType sp = sparseV.scalarProduct(bag[i]);
				if (std::norm(sp)>1e-8) return false;
			}
			return true;
		}

		void split(typename PsimagLite::Vector<SparseMatrixType>::Type& matrix,
		           const SparseMatrixType& matrix2) const
		{
			SizeType offset = 0;
			for (SizeType i=0;i<kspace_.size();i++) {
				SizeType blockSize = kspace_.blockSizes(i);
				if (blockSize==0) continue;
				std::cout<<"BLOCKSIZE="<<blockSize<<"\n";
				SparseMatrixType m(blockSize,blockSize);
				SizeType counter = 0;
				for (SizeType row=0;row<blockSize;row++) {
					m.setRow(row,counter);
					SizeType globalRow = row + offset;
					for (int k=matrix2.getRowPtr(globalRow);k<matrix2.getRowPtr(globalRow+1);k++) {
						ComplexOrRealType val = matrix2.getValue(k);
						if (std::norm(val)<1e-8) continue;
						SizeType globalCol = matrix2.getCol(k);
//						if (globalCol<offset) continue; // <<-- FIXME
						assert(globalCol>=offset);
						SizeType col = globalCol - offset;
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
				offset += blockSize;
			}
		}
//		{
//			SizeType counter = 0;
//			matrixA.resize(plusSector_,plusSector_);
//			for (SizeType i=0;i<plusSector_;i++) {
//				matrixA.setRow(i,counter);
//				for (int k=matrix.getRowPtr(i);k<matrix.getRowPtr(i+1);k++) {
//					SizeType col = matrix.getCol(k);
//					RealType val = matrix.getValue(k);
//					if (col<plusSector_) {
//						matrixA.pushCol(col);
//						matrixA.pushValue(val);
//						counter++;
//						continue;
//					}
//					if (!isAlmostZero(val)) {
//						PsimagLite::String s(__FILE__);
//						s += " Hamiltonian has no reflection symmetry.";
//						throw std::runtime_error(s.c_str());
//					}
//				}
//			}
//			matrixA.setRow(plusSector_,counter);

//			SizeType rank = matrix.row();
//			SizeType minusSector=rank-plusSector_;
//			matrixB.resize(minusSector,minusSector);
//			counter=0;
//			for (SizeType i=plusSector_;i<rank;i++) {
//				matrixB.setRow(i-plusSector_,counter);
//				for (int k=matrix.getRowPtr(i);k<matrix.getRowPtr(i+1);k++) {
//					SizeType col = matrix.getCol(k);
//					RealType val = matrix.getValue(k);
//					if (col>=plusSector_) {
//						matrixB.pushCol(col-plusSector_);
//						matrixB.pushValue(val);
//						counter++;
//						continue;
//					}
//					if (!isAlmostZero(val)) {
//						PsimagLite::String s(__FILE__);
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
		typename PsimagLite::Vector<SparseMatrixType>::Type matrixStored_;
		SizeType pointer_;
//		SparseMatrixType s_;
	}; // class TranslationSymmetry
} // namespace Dmrg

#endif  // TRANSLATION_SYMM_H
