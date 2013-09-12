#include "String.h"
PsimagLite::String license = "Copyright (c) 2009-2012, UT-Battelle, LLC\n"
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
#include "Concurrency.h"
#include "Engine.h"
#include "ProgramGlobals.h"

#include "Tj1Orb.h"
#include "Immm.h"
#include "HubbardOneOrbital.h"
#include "FeBasedSc.h"

#include "Geometry/Geometry.h"
#include "InternalProductStored.h"
#include "InputNg.h" // in PsimagLite
#include "ProgramGlobals.h"
#include "ContinuedFraction.h" // in PsimagLite 
#include "ContinuedFractionCollection.h" // in PsimagLite
#include "DefaultSymmetry.h"
#include "ReflectionSymmetry.h"
#include "TranslationSymmetry.h"
#include "Tokenizer.h"
#include "InputCheck.h"

using namespace LanczosPlusPlus;

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif
typedef std::complex<RealType> ComplexType;
typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::Geometry<RealType,ProgramGlobals> GeometryType;
typedef std::pair<SizeType,SizeType> PairType;
typedef PsimagLite::InputNg<InputCheck> InputNgType;

void fillOrbsOrSpin(PsimagLite::Vector<PairType>::Type& spinV,const PsimagLite::Vector<PsimagLite::String>::Type& strV)
{
	for (SizeType i=0;i<strV.size();i++) {
		PsimagLite::Vector<PsimagLite::String>::Type strV2;
		PsimagLite::tokenizer(strV[i],strV2,",");
		if (strV2.size()!=2)
			throw std::runtime_error("-o needs pairs\n");
		PairType spins;
		spins.first = atoi(strV2[0].c_str());
		spins.second = atoi(strV2[1].c_str());
		spinV.push_back(spins);
	}
}

template<typename ModelType>
SizeType maxOrbitals(const ModelType& model)
{
	SizeType res=0;
	for (SizeType i=0;i<model.geometry().numberOfSites();i++) {
		if (res<model.orbitals(i)) res=model.orbitals(i);
	}
	return res;
}

template<typename ModelType,typename SpecialSymmetryType>
void mainLoop2(ModelType& model,
			   InputNgType::Readable& io,
			   const GeometryType& geometry,
			   const PsimagLite::Vector<SizeType>::Type& gfV,
			   PsimagLite::Vector<SizeType>::Type& sites,
			   const PsimagLite::Vector<SizeType>::Type& cicjV,
               const PsimagLite::Vector<PairType>::Type& spins)
{
	typedef typename ModelType::BasisType BasisType;
	typedef Engine<ModelType,InternalProductStored,SpecialSymmetryType> EngineType;
	typedef typename EngineType::TridiagonalMatrixType TridiagonalMatrixType;

	EngineType engine(model,geometry.numberOfSites(),io);

	//! get the g.s.:
	RealType Eg = engine.gsEnergy();
	std::cout.precision(8);
	std::cout<<"Energy="<<Eg<<"\n";
	for (SizeType gfi=0;gfi<gfV.size();gfi++) {
		SizeType gf = gfV[gfi];
		io.read(sites,"TSPSites");
		if (sites.size()==0) throw std::runtime_error("No sites in input file!\n");
		if (sites.size()==1) sites.push_back(sites[0]);

		std::cout<<"#gf(i="<<sites[0]<<",j="<<sites[1]<<")\n";
		typedef PsimagLite::ContinuedFraction<TridiagonalMatrixType>
		ContinuedFractionType;
		typedef PsimagLite::ContinuedFractionCollection<ContinuedFractionType>
			ContinuedFractionCollectionType;

		PsimagLite::IoSimple::Out ioOut(std::cout);
		ContinuedFractionCollectionType cfCollection;
		SizeType norbitals = maxOrbitals(model);
		for (SizeType orb1=0;orb1<norbitals;orb1++) {
			for (SizeType orb2=0;orb2<norbitals;orb2++) {
				engine.spectralFunction(cfCollection,gf,sites[0],sites[1],spins,std::pair<SizeType,SizeType>(orb1,orb2));
			}
		}
		cfCollection.save(ioOut);
	}

	for (SizeType cicji=0;cicji<cicjV.size();cicji++) {
		SizeType cicj = cicjV[cicji];
		SizeType total = geometry.numberOfSites();
		PsimagLite::Matrix<typename SpecialSymmetryType::VectorType::value_type> cicjMatrix(total,total);
		SizeType norbitals = maxOrbitals(model);
		for (SizeType orb1=0;orb1<norbitals;orb1++) {
			for (SizeType orb2=0;orb2<norbitals;orb2++) {
				engine.twoPoint(cicjMatrix,cicj,spins,std::pair<SizeType,SizeType>(orb1,orb2));
				std::cout<<cicjMatrix;
			}
		}
	}
}

template<typename ModelType>
void mainLoop(InputNgType::Readable& io,
			  const GeometryType& geometry,
			  const PsimagLite::Vector<SizeType>::Type& gf,
			  PsimagLite::Vector<SizeType>::Type& sites,
			  const PsimagLite::Vector<SizeType>::Type& cicj,
              const PsimagLite::Vector<PairType>::Type& spins)
{
	typedef typename ModelType::ParametersModelType ParametersModelType;
	typedef typename ModelType::BasisType BasisType;

	// read model parameters
	ParametersModelType mp(io);

	std::cout<<mp;
	SizeType nup = 0;
	SizeType ndown = 0;

	try {
		io.read(nup,"TargetElectronsUp");
		io.read(ndown,"TargetElectronsDown");
	} catch (std::exception& e) {
		PsimagLite::Vector<RealType>::Type v;
		io.read(v,"TargetQuantumNumbers");
		if (v.size() < 2)
			throw PsimagLite::RuntimeError("TargetQuantumNumbers\n");
		nup = static_cast<SizeType>(v[0]*geometry.numberOfSites());
		ndown = static_cast<SizeType>(v[1]*geometry.numberOfSites());
	}

	//! Setup the Model
	ModelType model(nup,ndown,mp,geometry);

	int tmp = 0;
	try {
		io.readline(tmp,"UseTranslationSymmetry=");
	} catch(std::exception& e) {}

	bool useTranslationSymmetry = (tmp==1) ? true : false;

	try {
		io.readline(tmp,"UseReflectionSymmetry=");
	} catch(std::exception& e) {}

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
	PsimagLite::Vector<SizeType>::Type cicj,gf;
	PsimagLite::String file = "";
	PsimagLite::Vector<SizeType>::Type sites;
	PsimagLite::Vector<PairType>::Type spins(1,PairType(0,0));
	PsimagLite::Vector<PsimagLite::String>::Type str;
	InputCheck inputCheck;

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
			inputCheck.usage(argv[0]);
			return 1;
		}
	}

	if (file == "") {
		inputCheck.usage(argv[0]);
		return 1;
	}

	//! setup distributed parallelization
	SizeType npthreads = 1;
	ConcurrencyType concurrency(&argc,&argv,npthreads);

	//Setup the Geometry
	InputNgType::Writeable ioWriteable(file,inputCheck);
	InputNgType::Readable io(ioWriteable);
	GeometryType geometry(io);

	// print license
	if (ConcurrencyType::root()) std::cerr<<license;

	PsimagLite::String model("");
	io.readline(model,"Model=");

	if (model=="Tj1Orb") {
		mainLoop<Tj1Orb<RealType,GeometryType, InputNgType::Readable> >
		        (io,geometry,gf,sites,cicj,spins);
	} else if (model=="Immm") {
		mainLoop<Immm<RealType,GeometryType, InputNgType::Readable> >
		        (io,geometry,gf,sites,cicj,spins);
	} else if (model=="HubbardOneBand" ||
	           model=="HubbardOneBandExtended" ||
	           model=="SuperHubbardExtended") {
		mainLoop<HubbardOneOrbital<RealType,GeometryType, InputNgType::Readable> >
		        (io,geometry,gf,sites,cicj,spins);
	} else if (model=="FeAsBasedSc" || model=="FeAsBasedScExtended") {
		mainLoop<FeBasedSc<RealType,GeometryType, InputNgType::Readable> >
		        (io,geometry,gf,sites,cicj,spins);
	} else {
		std::cerr<<"No known model "<<model<<"\n";
		return 1;
	}
}

