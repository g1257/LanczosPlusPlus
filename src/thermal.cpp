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

	SizeType size() const { return eigs_.size(); }

	void multiplyRight(MatrixType& x,const MatrixType& a) const
	{
		SizeType n = a.n_row();
		SizeType m = a.n_col();
		assert(vecs_.n_row() == m);
		assert(vecs_.n_col() == m);
		assert(x.n_row() == n);
		assert(x.n_col() == m);
		assert(m > 0 && n > 0);
		psimag::BLAS::GEMM('N','N',n,m,m,1.0,&(a(0,0)),n,&(vecs_(0,0)),m,0.0,&(x(0,0)),n);
	}

	void multiplyLeft(MatrixType& x,const MatrixType& a) const
	{
		SizeType n = a.n_row();
		SizeType m = a.n_col();
		assert(vecs_.n_row() == n);
		psimag::BLAS::GEMM('C','N',n,m,n,1.0,&(vecs_(0,0)),n,&(a(0,0)),n,0.0,&(x(0,0)),n);
	}

	const RealType& eig(SizeType i) const
	{
		assert(i < eigs_.size());
		return eigs_[i];
	}

private:

	VectorSizeType sector_;
	VectorRealType eigs_;
	MatrixType vecs_;
};

typedef PsimagLite::Vector<OneSector*>::Type VectorOneSectorType;

struct ThermalOptions {
	ThermalOptions(PsimagLite::String operatorName_,
	               RealType beta_)
	    : operatorName(operatorName_),
	      beta(beta_)
	{}

	PsimagLite::String operatorName;
	RealType beta;
};

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

void findOperatorAndMatrix(MatrixType& a,
                           SizeType& jnd,
                           SizeType ind,
                           PsimagLite::String operatorName,
                           const VectorOneSectorType& sectors,
                           InputType& io)
{
	if (operatorName == "i") {
		jnd = ind;
		SizeType n = sectors[ind]->size();
		a.resize(n,n);
		a.setTo(0);
		for (SizeType i = 0; i < n; ++i)
			a(i,i) = 1.0;

		return;
	}

	VectorSizeType jndVector;
	io.read(jndVector,"#SectorDest");
	if (jndVector.size() == 0) return;
	jnd = findJnd(sectors,jndVector);

	io.readMatrix(a,"#Matrix");
}

RealType computeThisSector(SizeType ind,
                           const ThermalOptions& opt,
                           const VectorOneSectorType& sectors,
                           InputType& io)
{
	SizeType jnd = 0;
	MatrixType a;
	findOperatorAndMatrix(a,jnd,ind,opt.operatorName,sectors,io);

	SizeType n = a.n_row();
	if (n == 0) return 0.0;
	assert(n == sectors[ind]->size());

	SizeType m = a.n_col();
	assert(m == sectors[jnd]->size());
	if (m == 0) return 0.0;

	//Read operator 2   --> B
	MatrixType b = a;

	//Compute X^(s,s')_{n,n'} = \sum_{t,t'}U^{s*}_{n,t}A_{t,t'}^(s,s')U^{s'}_{t',n'}
	MatrixType x(n,m);
	computeX(x,a,*(sectors[ind]),*(sectors[jnd]));
	// Same for Y from B
	MatrixType y(n,m);
	computeX(y,b,*(sectors[ind]),*(sectors[jnd]));
	// Result is
	// \sum_{n,n'} X_{n,n'} Y_{n',n} exp(-i(E_n'-E_n)t))exp(-beta E_n)
	MatrixType z(n,n);
	psimag::BLAS::GEMM('N','C',n,n,m,1.0,&(x(0,0)),n,&(y(0,0)),n,0.0,&(z(0,0)),n);
	RealType sum = 0.0;
	for (SizeType i = 0; i < n; ++i) {
		sum += z(i,i) * exp(-opt.beta*sectors[ind]->eig(i));
	}

	std::cerr<<"--> "<<sum<<" n="<<n<<"\n";
	return sum;
}

void computeAverageFor(const ThermalOptions& opt,
                       const VectorOneSectorType& sectors,
                       InputType& io)
{
	RealType sum = 0.0;
	for (SizeType i = 0; i < sectors.size(); ++i)
		sum += computeThisSector(i,opt,sectors,io);

	std::cout<<"operator="<<opt.operatorName;
	std::cout<<" beta="<<opt.beta<<" sum="<<sum<<"\n";
}

void usage(char *name)
{
	std::cerr<<"USAGE: "<<name<<" -f file -c operator\n";
}

int main(int argc, char**argv)
{
	int opt = 0;
	PsimagLite::String operatorName;
	PsimagLite::String file;
	RealType beta = 0;
	while ((opt = getopt(argc, argv, "f:c:b:")) != -1) {
		switch (opt) {
		case 'c':
			operatorName = optarg;
			break;
		case 'f':
			file = optarg;
			break;
		case 'b':
			beta = atof(optarg);
			break;
		default: /* '?' */
			usage(argv[0]);
			return 1;
		}
	}

	if (file == "" || operatorName == "") {
		usage(argv[0]);
		return 1;
	}

	InputType io(file);
	SizeType total = 0;
	io.readline(total,"#TotalSectors=");
	VectorOneSectorType sectors(total);

	for (SizeType i = 0; i < sectors.size(); ++i) {
		sectors[i] = new OneSector(io);
		sectors[i]->info(std::cout);
	}

	ThermalOptions options(operatorName,beta);
	io.rewind();
	computeAverageFor(options,sectors,io);

	for (SizeType i = 0; i < sectors.size(); ++i)
		delete sectors[i];
}

