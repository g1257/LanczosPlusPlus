/*
Copyright (c) 2014-2016, UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/

#ifndef LANCZOS_REDUCED_DM_H
#define LANCZOS_REDUCED_DM_H
#include "Vector.h"

namespace LanczosPlusPlus {

//template<typename GeometryType>
class ReducedDensityMatrix {

	typedef double RealType;

public:

	ReducedDensityMatrix()
	{
		build();
		w_ = rdm_;
		diag(diag_,eigs_,'V');
	}

	void printAll(std::ostream& os) const
	{
		os<<"Reduced Density Matrix\n";
		os<<rdm_;
		os<<"Eigenvectors of Reduced Density Matrix\n";
		os<<w_;
		os<<"Eigenvalues of Reduced Density Matrix\n";
		os<<eigs_;
	}

private:

	void build()
	{
		for (SizeType alpha = 0; alpha < ns; ++alpha) {
			for (SizeType alphaP = 0; alphaP < ns; ++alphaP) {
				rdm_(alpha,alphaP) = rdm(alpha,alphaP);
			}
		}
	}

	RealType rdm(SizeType alpha, SizeType alphaP) const
	{
		ComplexOrRealType sum = 0.0;
		for (SizeType beta = 0; beta < ne; ++beta) {
			// alpha + beta --> i
			SizeType i = pack(alpha,beta);
			// alphaP + beta --> j
			SizeType j = pack(alphaP,beta);
			sum += std::conj(psi_[i])*psi_[j];
		}

		return sum;
	}

	MatrixType rdm_;
	MatrixType w_;
	VectorRealType eigs_;
}; // class ReducedDensityMatrix

} // namespace LanczosPlusPlus
#endif

