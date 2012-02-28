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

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" [-g -i i -j j] -f filename\n";
}

RealType proc(const SparseMatrixType& hti,const SparseMatrixType& htd,const RealType& time)
{
	PsimagLite::Matrix<RealType> fhti;
	crsMatrixToFullMatrix(fhti,hti);
	std::vector<RealType> eigs(fhti.n_row());
	std::cerr<<"Diag. size="<<eigs.size()<<"\n";
	diag(fhti,eigs,'V');
	std::cout<<"Lowest energy: time independent="<<eigs[0]<<" "<<(2.0*eigs[0])<<"\n";


	PsimagLite::Matrix<RealType> fhtd;
	crsMatrixToFullMatrix(fhtd,htd);

	diag(fhtd,eigs,'V');
	std::cout<<"Lowest energy: time dependent="<<eigs[0]<<" "<<(2.0*eigs[0])<<" time="<<time<<"\n";

	RealType sum2 = 0;

	size_t n = fhtd.n_row();
	for (size_t i=0;i<n;i++) {
		RealType sum = 0.0;
		for (size_t k=0;k<n;k++) {
			for (size_t k2=0;k2<n;k2++) {
				sum += std::conj(fhtd(k,i))*fhtd(k2,i)*std::conj(fhti(k,0))*fhti(k2,0);
			}
		}
		sum2 += eigs[i] * sum;
	}
	std::cout<<"<psi(0)|S S lambda_i | psi(0)>="<<sum2<<" "<<(2.0*sum2)<<" time="<<time<<"\n";
	return (2.0*sum2);
}

int main(int argc,char *argv[])
{
	int opt = 0;
	bool gf = false;
	std::string file = "";
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


	std::vector<RealType> qns;
	io.read(qns,"TargetQuantumNumbers");
	if (qns.size()<2) throw std::runtime_error("HubbardLanczos::ctor(...)\n");
	size_t nup=size_t(geometry.numberOfSites()*qns[0]);
	size_t ndown=size_t(geometry.numberOfSites()*qns[1]);

	ParametersModelType mpTimeIndepedent = mp;
	mpTimeIndepedent.omegaTime = 0.0;

	ModelType modelTimeIndependent(nup,ndown,mpTimeIndepedent,geometry);
	SparseMatrixType hTimeIndependent;
	modelTimeIndependent.setupHamiltonian(hTimeIndependent);

	//! Setup the Models
	size_t numberOfTimes = 11;
	RealType deltaTime = 0.1;
	RealType omega = 0.8;
	std::vector<RealType> v(numberOfTimes);
	for (size_t i=0;i<numberOfTimes;i++) {
		RealType time = i*deltaTime;
		mp.omegaTime = omega * time;

		ModelType modelTimeDependent(nup,ndown,mp,geometry);
		SparseMatrixType hTimeDependent;
		modelTimeDependent.setupHamiltonian(hTimeDependent);
		v[i] = proc(hTimeIndependent,hTimeDependent,time);
	}
	for (size_t i=0;i<v.size();i++) {
		RealType time = i*deltaTime;
		std::cout<<time<<" "<<v[i]<<"\n";
	}
}

