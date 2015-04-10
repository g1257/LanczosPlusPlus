#include "Vector.h"
#include "IoSimple.h"
#include "Matrix.h"
#include "BLAS.h"

typedef double RealType;
typedef PsimagLite::Matrix<RealType> MatrixType;
typedef PsimagLite::Vector<RealType>::Type VectorRealType;
typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
typedef PsimagLite::IoSimple::In InputType;

class OneSector {

public:

	OneSector(InputType& io)
	{
		io.read(sector_,"#SectorSource");
		io.read(eigs_,"#Eigenvalues");
		io.readMatrix(vecs_,"#Eigenvectors");
	}

	bool isSector(const VectorSizeType& jndVector) const
	{
		return (jndVector == sector_);
	}

	void info(std::ostream& os) const
	{
		os<<"sector\n";
		os<<sector_;
		os<<"eigs.size()="<<eigs_.size()<<"\n";
		os<<"vecs="<<vecs_.n_row()<<"x"<<vecs_.n_col()<<"\n";
	}

	void multiplyRight(MatrixType& x,const MatrixType& a) const
	{
		SizeType n = a.n_row();
		psimag::BLAS::GEMM('N','N',n,n,n,1.0,&(a(0,0)),n,&(vecs_(0,0)),n,0.0,&(x(0,0)),n);
	}

	void multiplyLeft(MatrixType& x,const MatrixType& a) const
	{
		SizeType n = a.n_row();
		psimag::BLAS::GEMM('N','N',n,n,n,1.0,&(vecs_(0,0)),n,&(a(0,0)),n,0.0,&(x(0,0)),n);
	}

private:

	VectorSizeType sector_;
	VectorRealType eigs_;
	MatrixType vecs_;
};

typedef PsimagLite::Vector<OneSector*>::Type VectorOneSectorType;

//Compute X^(s,s')_{n,n'} = \sum_{t,t'}U^{s*}_{n,t}A_{t,t'}^(s,s')U^{s'}_{t',n'}
void computeX(MatrixType& x,
              const MatrixType& a,
              const OneSector& sectorSrc,
              const OneSector& sectorDest)
{
	MatrixType tmp(x.n_row(),x.n_col());
	sectorDest.multiplyRight(tmp,a);
	sectorSrc.multiplyLeft(x,tmp);
}

SizeType findJnd(const VectorOneSectorType& sectors, const VectorSizeType& jndVector)
{
	for (SizeType i = 0; i < sectors.size(); ++i) {
		if (sectors[i]->isSector(jndVector)) return i;
	}

	PsimagLite::String str(__FILE__);
	str += " " + ttos(__LINE__) + "\n";
	str += "findJnd can't find sector index\n";
	throw PsimagLite::RuntimeError(str);
}

void computeThisSector(SizeType ind, const VectorOneSectorType& sectors, InputType& io)
{
	VectorSizeType jndVector;
	io.read(jndVector,"#SectorDest");
	if (jndVector.size() == 0) return;
	SizeType jnd = findJnd(sectors,jndVector);

	MatrixType a;
	io.readMatrix(a,"#Matrix");
	SizeType n = a.n_row();
	if (n == 0) return;
	//Read operator 2   --> B
	MatrixType b = a;

	//Compute X^(s,s')_{n,n'} = \sum_{t,t'}U^{s*}_{n,t}A_{t,t'}^(s,s')U^{s'}_{t',n'}
	MatrixType x(n,n);
	computeX(x,a,*(sectors[ind]),*(sectors[jnd]));
	// Same for Y from B
	MatrixType y(n,n);
	computeX(y,b,*(sectors[ind]),*(sectors[jnd]));
	// Result is
	// \sum_{n,n'} X_{n,n'} Y_{n',n} exp(-i(E_n'-E_n)t))exp(-beta E_n)
}

void computeAverageFor(const VectorOneSectorType& sectors, InputType& io)
{
	for (SizeType i = 0; i < sectors.size(); ++i)
		computeThisSector(i,sectors,io);
}

int main(int argc, char**argv)
{
	if (argc < 2) return 1;

	InputType io(argv[1]);
	SizeType total = 0;
	io.readline(total,"#TotalSectors=");
	VectorOneSectorType sectors(total);

	for (SizeType i = 0; i < sectors.size(); ++i) {
		sectors[i] = new OneSector(io);
		sectors[i]->info(std::cout);
	}

	InputType io2(argv[1]);
	computeAverageFor(sectors,io2);

	for (SizeType i = 0; i < sectors.size(); ++i)
		delete sectors[i];
}

