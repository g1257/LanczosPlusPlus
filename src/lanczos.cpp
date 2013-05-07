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
#include "DefaultSymmetry.h"
#include "ReflectionSymmetry.h"
#include "TranslationSymmetry.h"
#include "Tokenizer.h"

using namespace LanczosPlusPlus;

typedef double RealType;
typedef std::complex<RealType> ComplexType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::Geometry<RealType,ProgramGlobals> GeometryType;
typedef PsimagLite::IoSimple::In IoInputType;
typedef std::pair<size_t,size_t> PairType;

void fillOrbsOrSpin(std::vector<PairType>& spinV,const std::vector<std::string>& strV)
{
	for (size_t i=0;i<strV.size();i++) {
		std::vector<std::string> strV2;
		PsimagLite::tokenizer(strV[i],strV2,",");
		if (strV2.size()!=2)
			throw std::runtime_error("-o needs pairs\n");
		PairType spins;
		spins.first = atoi(strV2[0].c_str());
		spins.second = atoi(strV2[1].c_str());
		spinV.push_back(spins);
	}
}

std::pair<size_t,size_t> readElectrons(PsimagLite::IoSimple::In& io,size_t nsites)
{
	int nup = -1;
	int ndown = -1;
	try {
		io.readline(nup,"TargetElectronsUp=");
		io.readline(ndown,"TargetElectronsDown=");
	} catch (std::exception& e)
	{
		nup = ndown = -1;
		io.rewind();
	}

	std::vector<RealType> v;
	try {
		io.read(v,"TargetQuantumNumbers");
	} catch (std::exception& e)
	{
		v.resize(0);
		io.rewind();
	}

	if (nup<0 && v.size()==0) {
		std::string str("Either TargetElectronsUp/Down or TargetQuantumNumbers is need\n");
		throw std::runtime_error(str.c_str());
	}

	if (nup>=0 && v.size()>0) {
		std::string str("Having both TargetElectronsUp/Down and TargetQuantumNumbers is an error\n");
		throw std::runtime_error(str.c_str());
	}

	if (nup>=0) return std::pair<size_t,size_t>(nup,ndown);

	if (v.size()<2) {
		std::string str("Incorrect TargetQuantumNumbers line\n");
		throw std::runtime_error(str.c_str());
	}
	nup = size_t(v[0]*nsites);
	ndown = size_t(v[1]*nsites);
	return std::pair<size_t,size_t>(nup,ndown);
}

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" [-g -c] -f filename\n";
}

template<typename ModelType>
size_t maxOrbitals(const ModelType& model)
{
	size_t res=0;
	for (size_t i=0;i<model.geometry().numberOfSites();i++) {
		if (res<model.orbitals(i)) res=model.orbitals(i);
	}
	return res;
}

template<typename ModelType,typename SpecialSymmetryType>
void mainLoop2(ModelType& model,
			   IoInputType& io,
			   const GeometryType& geometry,
			   const std::vector<size_t>& gfV,
			   std::vector<size_t>& sites,
			   const std::vector<size_t>& cicjV,
               const std::vector<PairType>& spins)
{
	typedef typename ModelType::BasisType BasisType;
	typedef Engine<ModelType,InternalProductStored,SpecialSymmetryType,ConcurrencyType> EngineType;
	typedef typename EngineType::TridiagonalMatrixType TridiagonalMatrixType;

	EngineType engine(model,geometry.numberOfSites(),io);

	//! get the g.s.:
	RealType Eg = engine.gsEnergy();
	std::cout.precision(8);
	std::cout<<"Energy="<<Eg<<"\n";
	for (size_t gfi=0;gfi<gfV.size();gfi++) {
		size_t gf = gfV[gfi];
		io.read(sites,"TSPSites");
		if (sites.size()==0) throw std::runtime_error("No sites in input file!\n");
		if (sites.size()==1) sites.push_back(sites[0]);

		std::cout<<"#gf(i="<<sites[0]<<",j="<<sites[1]<<")\n";
		typedef PsimagLite::ContinuedFraction<RealType,TridiagonalMatrixType>
		ContinuedFractionType;
		typedef PsimagLite::ContinuedFractionCollection<ContinuedFractionType>
			ContinuedFractionCollectionType;

		PsimagLite::IoSimple::Out ioOut(std::cout);
		ContinuedFractionCollectionType cfCollection;
		size_t norbitals = maxOrbitals(model);
		for (size_t orb1=0;orb1<norbitals;orb1++) {
			for (size_t orb2=0;orb2<norbitals;orb2++) {
				engine.spectralFunction(cfCollection,gf,sites[0],sites[1],spins,std::pair<size_t,size_t>(orb1,orb2));
			}
		}
		cfCollection.save(ioOut);
	}

	for (size_t cicji=0;cicji<cicjV.size();cicji++) {
		size_t cicj = cicjV[cicji];
		size_t total = geometry.numberOfSites();
		PsimagLite::Matrix<typename SpecialSymmetryType::VectorType::value_type> cicjMatrix(total,total);
		size_t norbitals = maxOrbitals(model);
		for (size_t orb1=0;orb1<norbitals;orb1++) {
			for (size_t orb2=0;orb2<norbitals;orb2++) {
				engine.twoPoint(cicjMatrix,cicj,spins,std::pair<size_t,size_t>(orb1,orb2));
				std::cout<<cicjMatrix;
			}
		}
	}
}

template<typename ModelType>
void mainLoop(IoInputType& io,
			  const GeometryType& geometry,
			  const std::vector<size_t>& gf,
			  std::vector<size_t>& sites,
			  const std::vector<size_t>& cicj,
              const std::vector<PairType>& spins)
{
	typedef typename ModelType::ParametersModelType ParametersModelType;
	typedef typename ModelType::BasisType BasisType;

	// read model parameters
	ParametersModelType mp(io);

	std::cout<<mp;
	std::pair<size_t,size_t> nupndown = readElectrons(io,geometry.numberOfSites());

	//! Setup the Model
	ModelType model(nupndown.first,nupndown.second,mp,geometry);

	int tmp = 0;
	try {
		io.readline(tmp,"UseTranslationSymmetry=");
	} catch(std::exception& e) {}
	io.rewind();
	bool useTranslationSymmetry = (tmp==1) ? true : false;

	try {
		io.readline(tmp,"UseReflectionSymmetry=");
	} catch(std::exception& e) {}
	io.rewind();
	bool useReflectionSymmetry = (tmp==1) ? true : false;

	if (useTranslationSymmetry) {
		mainLoop2<ModelType,TranslationSymmetry<GeometryType,BasisType> >(model,io,geometry,gf,sites,cicj,spins);
	} else if (useReflectionSymmetry) {
		mainLoop2<ModelType,ReflectionSymmetry<GeometryType,BasisType> >(model,io,geometry,gf,sites,cicj,spins);
	} else {
		mainLoop2<ModelType,DefaultSymmetry<GeometryType,BasisType> >(model,io,geometry,gf,sites,cicj,spins);
	}
}


int main(int argc,char *argv[])
{
	int opt = 0;
	std::vector<size_t> cicj,gf;
	std::string file = "";
	std::vector<size_t> sites;
	std::vector<PairType> spins(1,PairType(0,0));
	std::vector<std::string> str;

	while ((opt = getopt(argc, argv, "g:c:f:s:")) != -1) {
		switch (opt) {
		case 'g':
			gf.push_back(ProgramGlobals::operator2id(optarg));
			break;
		case 'f':
			file = optarg;
			break;
		case 'c':
			cicj.push_back(ProgramGlobals::operator2id(optarg));
			break;
		case 's':
			PsimagLite::tokenizer(optarg,str,";");
			fillOrbsOrSpin(spins,str);
			str.clear();
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
		mainLoop<Tj1Orb<RealType,GeometryType> >(io,geometry,gf,sites,cicj,spins);
	} else if (model=="Immm") {
		mainLoop<Immm<RealType,GeometryType> >(io,geometry,gf,sites,cicj,spins);
	} else if (model=="HubbardOneBand") {
		mainLoop<HubbardOneOrbital<RealType,GeometryType> >(io,geometry,gf,sites,cicj,spins);
	} else if (model=="FeAsBasedSc") {
		mainLoop<FeBasedSc<RealType,GeometryType> >(io,geometry,gf,sites,cicj,spins);
	} else {
		std::cerr<<"No known model "<<model<<"\n";
		return 1;
	}
}

