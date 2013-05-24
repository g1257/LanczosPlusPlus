/* DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
 * This driver program was written by configure.pl
 * Lanczos++ (v1.0) by G.A.*/
#include <unistd.h>
#include <cstdlib>
#include <getopt.h>
#include "ConcurrencySerial.h"
#include "Engine.h"
#include "HubbardOneOrbital.h"
#include "ParametersModelHubbard.h"
#include "Geometry.h"
#include "InternalProductStored.h"
#include "IoSimple.h" // in PsimagLite
#include "ProgramGlobals.h"

#include "ReflectionSymmetry.h"

using namespace LanczosPlusPlus;

typedef double RealType;
typedef std::complex<RealType> ComplexType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::Geometry<RealType,ProgramGlobals> GeometryType;
typedef ParametersModelHubbard<RealType> ParametersModelType;
typedef PsimagLite::IoSimple::In IoInputType;
typedef HubbardOneOrbital<RealType,ParametersModelType,GeometryType> ModelType;
typedef typename ModelType::BasisType BasisType;
typedef ReflectionSymmetry<GeometryType,BasisType> ReflectionSymmetryType;
typedef InternalProductStored<ModelType,ReflectionSymmetryType> InternalProductType;
typedef Engine<ModelType,InternalProductType,ConcurrencyType> EngineType;
typedef typename InternalProductType::SparseMatrixType SparseMatrixType;
typedef std::complex<RealType> ComplexType;

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" [-g -i i -j j] -f filename\n";
}

ComplexType proc(PsimagLite::Vector<ComplexType>::Type& psi,
		 const SparseMatrixType& hmat,
		 const SparseMatrixType& hIntegral,
		 const RealType& tau)
{
	PsimagLite::Matrix<ComplexType> fhmat;
	crsMatrixToFullMatrix(fhmat,hmat);

	SizeType n = fhmat.n_row();
	PsimagLite::Matrix<ComplexType> fhIntegral(n,n);
	crsMatrixToFullMatrix(fhIntegral,hIntegral);
	ComplexType val= ComplexType(0,-tau);
	fhIntegral = (val * fhIntegral);

	exp(fhIntegral);

	PsimagLite::Vector<ComplexType>::Type psiNew(n);
	for (SizeType i=0;i<n;i++) {
		ComplexType sum = 0;
		for (SizeType j=0;j<n;j++)  sum += fhIntegral(i,j) * psi[j];
		psiNew[i] = sum;
	}

	PsimagLite::Vector<ComplexType>::Type v = fhmat * psiNew;

	//std::cout<<"<psi(0)|S S lambda_i | psi(0)>="<<sum2<<" "<<(2.0*sum2)<<" time="<<time<<"\n";
	psi = psiNew;
	return psiNew * v;
}

int main(int argc,char *argv[])
{
	int opt = 0;
	bool gf = false;
	PsimagLite::String file = "";
	while ((opt = getopt(argc, argv, "gf:")) != -1) {
		switch (opt) {
		case 'g':
			gf = true;
			break;
		case 'f':
			file = optarg;
			break;
		default: /* '?' */
			usage(argv[0]);
			return 1;
		}
	}
	if (file == "") {
		usage(argv[0]);
		return 1;
	}
	//! setup distributed parallelization
	ConcurrencyType concurrency(argc,argv);

	//Setup the Geometry
	IoInputType io(file);
	GeometryType geometry(io);

	// read model parameters
	ParametersModelType mp(io);


	PsimagLite::Vector<RealType>::Type qns;
	io.read(qns,"TargetQuantumNumbers");
	if (qns.size()<2) throw std::runtime_error("HubbardLanczos::ctor(...)\n");
	SizeType nup=SizeType(geometry.numberOfSites()*qns[0]);
	SizeType ndown=SizeType(geometry.numberOfSites()*qns[1]);

	ParametersModelType mpTimeIndepedent = mp;

	mpTimeIndepedent.timeFactor = 1.0;

	ModelType modelTimeIndependent2(nup,ndown,mpTimeIndepedent,geometry);
	SparseMatrixType psiM;
	modelTimeIndependent2.setupHamiltonian(psiM);

	PsimagLite::Matrix<RealType> fpsi;
	crsMatrixToFullMatrix(fpsi,psiM);
	PsimagLite::Vector<RealType>::Type eigs2(fpsi.n_row());
	diag(fpsi,eigs2,'V');
	PsimagLite::Vector<ComplexType>::Type psi(eigs2.size());
	for (SizeType i=0;i<psi.size();i++) psi[i] = fpsi(i,0);

	//! Setup the Models
	SizeType numberOfTimes = 510;
	RealType deltaTime = 0.01;
	RealType omega = 0.8;
	PsimagLite::Vector<ComplexType>::Type v(numberOfTimes);
	for (SizeType i=0;i<numberOfTimes;i++) {
		RealType time = i*deltaTime;

		mp.timeFactor = cos(omega*time);
		ModelType modelH(nup,ndown,mp,geometry);
		SparseMatrixType hmat;
		modelH.setupHamiltonian(hmat);

//		mp.timeFactor = 0;
//		if (time>0) mp.timeFactor = sin(omega*time)/(omega*deltaTime);
//		ModelType modelHIntegral(nup,ndown,mp,geometry);
//		SparseMatrixType hIntegral;
//		modelHIntegral.setupHamiltonian(hIntegral);
		v[i] = proc(psi,hmat,hmat,deltaTime);
		std::cerr<<time<<" "<<v[i]<<"\n";
	}

	for (SizeType i=0;i<v.size();i++) {
		RealType time = i*deltaTime;
		std::cout<<time<<" "<<std::real(v[i])<<"\n";
	}
}

