#include <iostream>
#include <cmath>
#include "../../PsimagLite/src/Matrix.h"
#include <cassert>
#include <vector>

template<typename ComplexOrRealType>
class Hamiltonian {

public:

	typedef std::vector<ComplexOrRealType> VectorType;

	Hamiltonian(int statesS, int statesL, int nsites)
		: statesS_(statesS), 
		statesL_(statesL), 
		nsites_(nsites),
		sPerSite_(nsites),
		lPerSite_(nsites)
	{}

	void fill()
	{
		VectorType rowVector(statesS_*statesL_);
		for (int idS = 0; idS < statesS_; ++idS) {
			for (int idL = 0; idL < statesL_; ++idL) {
				fillMatrixRow(rowVector, idS, idL);
			}
		}
	}

private:

	int fillMatrixRow(VectorType& rowVector, int idS, int idL)
	{
		indexToVector(sPerSite_, idS);
		indexToVector(lPerSite_, idL);
		for (int i = 0; i < nsites_; ++i) {
			for (int j = i + 1; j < nsites_; ++j) {
				for (int which0 = 0; which0 < 3; ++which0) { // +- -+ or zz
					for (int which1 = 0; which1 < 3; ++which1) { // +- -+ or zz
						ComplexOrRealType valueS = 0;
						int braS = createOneTerm(valueS, idS, i, j, which0, 0);
						if (braS < 0) continue;
						
						ComplexOrRealType valueL = 0;
						int braL = createOneTerm(valueL, idL, i, j, which1, 1);
						if (braL < 0) continue;
						int col = packSandL(braS, braL);
						rowVector[col] += valueS*valueL;
					}
				}
			}
		}
	}

	int createOneTerm(ComplexOrRealType& value, int id, int i, int j, int which, int what)
	{
		std::vector<int>& v = (what == 0) ? sPerSite_ : lPerSite_;
		int mi = v[i];
		int mj = v[j];
		switch (which) {
		case 0:
			value = 1; // = 0.5*sqrt(2)^2
			return createSpSm(i, v[i], j, v[j]);
			break;
		case 1:
			value = 1; // = 0.5*sqrt(2)^2
			return createSmSp(i, v[i], j, v[j]);
			break;
		case 2:
			value = (v[i] - 1)*(v[j] - 1);
			return id;
		default:
			throw std::runtime_error("createOneTerm\n");
		}
	}

	int createSpSm(int i, int mi, int j, int mj)
	{
		if (mi == 2) return -1;
		if (mj == 0) return -1;
		return createOneState(i, mi + 1, j, mj - 1);
	}

	int createSmSp(int id, int i, int mi, int j, int mj)
	{
		if (mj == 2) return -1;
		if (mi == 0) return -1;
		return createOneState(i, mi - 1, j, mj + 1);
	}

	int createOneState(int ind, int mind, int jnd, int mjnd)
	{
		int savedInd = sPerSite_[ind];
		sPerSite_[ind] = mind;
		int savedJnd = sPerSite_[jnd];
		sPerSite_[jnd] = mjnd;
		int braIndex = vectorToIndex(sPerSite_);
		sPerSite_[ind] = savedInd;
		sPerSite_[jnd] = savedJnd;
		return braIndex;
	}
	
	int packSandL(int idS, int idL) const
	{
		assert(idS < statesS_);
		return idS + idL*statesS_;
	}

	int statesS_;
	int statesL_;
	int nsites_;
	std::vector<int> sPerSite_;
	std::vector<int> lPerSite_;
};

int main(int argc, char* argv[])
{
	if (argc != 2) {
		std::cerr<<"USAGE: "<<argv[0]<<" nsites\n";
		return 1;
	}

	int nsites = atoi(argv[1]);
	int statesS = pow(3, nsites);
	int statesL = pow(3, nsites);
	int states = statesS*statesL;
	PsimagLite::Matrix<double> matrix(states, states);

	Hamiltonian<double> h(statesS, statesL, nsites);

	h.fill();
}

