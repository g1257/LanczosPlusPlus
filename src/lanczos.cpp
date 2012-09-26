#include <string>
std::string license = "Copyright (c) 2009-2012, UT-Battelle, LLC\n"
"All rights reserved\n"
"\n"
"[Lanczos++, Version 1.0.0]\n"
"\n"
"*********************************************************\n"
"THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND\n"
"CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED\n"
"WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
"WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n"
"PARTICULAR PURPOSE ARE DISCLAIMED. \n"
"\n"
"Please see full open source license included in file LICENSE.\n"
"*********************************************************\n"
"\n";

#include <unistd.h>
#include <cstdlib>
#include <getopt.h>
#include "ConcurrencySerial.h"
#include "Engine.h"
#include "ProgramGlobals.h"

#include "Tj1Orb.h"
#include "Immm.h"
#include "HubbardOneOrbital.h"
#include "FeBasedSc.h"

#include "Geometry.h"
#include "InternalProductStored.h"
#include "IoSimple.h" // in PsimagLite
#include "ProgramGlobals.h"
#include "ContinuedFraction.h" // in PsimagLite 
#include "ContinuedFractionCollection.h" // in PsimagLite
#include "ReflectionSymmetry.h"
#include "Split.h"

using namespace LanczosPlusPlus;

typedef double RealType;
typedef std::complex<RealType> ComplexType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::Geometry<RealType,ProgramGlobals> GeometryType;
typedef PsimagLite::IoSimple::In IoInputType;

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" [-g -c] -f filename\n";
}


template<typename ModelType>
void mainLoop(IoInputType& io,const GeometryType& geometry,size_t gf,std::vector<size_t>& sites,size_t cicj)
{
	typedef typename ModelType::ParametersModelType ParametersModelType;
	typedef typename ModelType::BasisType BasisType;
	typedef ReflectionSymmetry<GeometryType,BasisType> ReflectionSymmetryType;
	typedef InternalProductStored<ModelType,ReflectionSymmetryType> InternalProductType;
	typedef Engine<ModelType,InternalProductType,ConcurrencyType> EngineType;
	typedef typename EngineType::TridiagonalMatrixType TridiagonalMatrixType;

	// read model parameters
	ParametersModelType mp(io);

	size_t nup = 0;
	size_t ndown = 0;
	io.readline(nup,"TargetElectronsUp=");
	io.readline(ndown,"TargetElectronsDown=");

	//! Setup the Model
	ModelType model(nup,ndown,mp,geometry);

	EngineType engine(model,geometry.numberOfSites(),io);

	//! get the g.s.:
	RealType Eg = engine.gsEnergy();
	std::cout.precision(8);
	std::cout<<"Energy="<<Eg<<"\n";
	if (gf!=ProgramGlobals::OPERATOR_NIL) {
		io.read(sites,"TSPSites");
		if (sites.size()==0) throw std::runtime_error("No sites in input file!\n");
		if (sites.size()==1) sites.push_back(sites[0]);

		std::cout<<"#gf(i="<<sites[0]<<",j="<<sites[1]<<")\n";
		typedef PsimagLite::ContinuedFraction<RealType,TridiagonalMatrixType>
		ContinuedFractionType;
		typedef PsimagLite::ContinuedFractionCollection<ContinuedFractionType>
			ContinuedFractionCollectionType;

		ContinuedFractionCollectionType cfCollection;
		engine.spectralFunction(cfCollection,gf,sites[0],sites[1],ModelType::SPIN_UP,std::pair<size_t,size_t>(0,0));

		PsimagLite::IoSimple::Out ioOut(std::cout);
		cfCollection.save(ioOut);
	}

	if (cicj!=ProgramGlobals::OPERATOR_NIL) {
		size_t total = geometry.numberOfSites();
		PsimagLite::Matrix<RealType> cicjMatrix(total,total);
		size_t norbitals = model.orbitals();
		for (size_t i=0;i<norbitals;i++) {
			engine.twoPoint(cicjMatrix,cicj,ModelType::SPIN_UP,std::pair<size_t,size_t>(i,i));
			std::cout<<cicjMatrix;
		}
	}
}

int main(int argc,char *argv[])
{
	int opt = 0;
	size_t gf = ProgramGlobals::OPERATOR_NIL;
	std::string file = "";
	std::vector<size_t> sites;
	size_t cicj=ProgramGlobals::OPERATOR_NIL;
	while ((opt = getopt(argc, argv, "g:c:f:")) != -1) {
		switch (opt) {
		case 'g':
			gf = ProgramGlobals::operator2id(optarg);
			break;
		case 'f':
			file = optarg;
			break;
		case 'c':
			cicj = ProgramGlobals::operator2id(optarg);
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

	// print license
	if (concurrency.root()) std::cerr<<license;

	std::string model("");
	io.readline(model,"Model=");

	if (model=="Tj1Orb") {
		mainLoop<Tj1Orb<RealType,GeometryType> >(io,geometry,gf,sites,cicj);
	} else if (model=="Immm") {
		mainLoop<Immm<RealType,GeometryType> >(io,geometry,gf,sites,cicj);
	} else if (model=="HubbardOneBand") {
		mainLoop<HubbardOneOrbital<RealType,GeometryType> >(io,geometry,gf,sites,cicj);
	} else if (model=="FeBasedSc") {
		mainLoop<FeBasedSc<RealType,GeometryType> >(io,geometry,gf,sites,cicj);
	}
}

