#include "AllocatorCpu.h"
#include "Version.h"
#include "../../PsimagLite/src/Version.h"
PsimagLite::String license = "Copyright (c) 2009-2012, UT-Battelle, LLC\n"
                             "All rights reserved\n"
                             "\n"
                             "[Lanczos++, Version 1.0]\n"
                             "\n"
                             "-------------------------------------------------------------\n"
                             "THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND\n"
                             "CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED\n"
                             "WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
                             "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n"
                             "PARTICULAR PURPOSE ARE DISCLAIMED. \n"
                             "\n"
                             "Please see full open source license included in file LICENSE.\n"
                             "-------------------------------------------------------------\n"
                             "\n";

#include <unistd.h>
#include <cstdlib>
#include <getopt.h>
#define USE_PTHREADS_OR_NOT_NG
#include "Concurrency.h"
#include "Engine.h"
#include "ProgramGlobals.h"
#include "ModelSelector.h"
#include "Geometry/Geometry.h"
#include "InternalProductOnTheFly.h"
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
#include "ReducedDensityMatrix.h"

using namespace LanczosPlusPlus;

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif

typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::InputNg<InputCheck> InputNgType;
typedef std::pair<SizeType,SizeType> PairSizeType;

struct LanczosOptions {

	LanczosOptions()
	    : split(-1),spins(1,PairSizeType(0,0)),extendedStatic("")
	{}

	int split;
	PsimagLite::Vector<SizeType>::Type cicj;
	PsimagLite::Vector<SizeType>::Type gf;
	PsimagLite::Vector<SizeType>::Type sites;
	PsimagLite::Vector<PairSizeType>::Type spins;
	PsimagLite::String extendedStatic;

}; // struct LanczosOptions

void fillOrbsOrSpin(PsimagLite::Vector<PairSizeType>::Type& spinV,
                    const PsimagLite::Vector<PsimagLite::String>::Type& strV)
{
	for (SizeType i=0;i<strV.size();i++) {
		PsimagLite::Vector<PsimagLite::String>::Type strV2;
		PsimagLite::tokenizer(strV[i],strV2,",");
		if (strV2.size()!=2)
			throw std::runtime_error("-o needs pairs\n");
		PairSizeType spins;
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

template<typename ModelType,
         typename SpecialSymmetryType,
         template<typename,typename> class InternalProductTemplate>
void mainLoop3(const ModelType& model,
               InputNgType::Readable& io,
               LanczosOptions& lanczosOptions)
{
	typedef typename ModelType::GeometryType GeometryType;
	typedef typename GeometryType::ComplexOrRealType ComplexOrRealType;
	typedef Engine<ModelType,InternalProductTemplate,SpecialSymmetryType> EngineType;
	typedef typename EngineType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef typename EngineType::VectorSizeType VectorSizeType;

	const GeometryType& geometry = model.geometry();
	EngineType engine(model,geometry.numberOfSites(),io);

	//! get the g.s.:
	RealType Eg = engine.gsEnergy();
	std::cout.precision(8);
	std::cout<<"Energy="<<Eg<<"\n";
	for (SizeType gfi=0;gfi<lanczosOptions.gf.size();gfi++) {
		SizeType gfI = lanczosOptions.gf[gfi];
		io.read(lanczosOptions.sites,"TSPSites");
		if (lanczosOptions.sites.size()==0)
			throw std::runtime_error("No sites in input file!\n");
		if (lanczosOptions.sites.size()==1)
			lanczosOptions.sites.push_back(lanczosOptions.sites[0]);

		std::cout<<"#gf(i="<<lanczosOptions.sites[0]<<",j=";
		std::cout<<lanczosOptions.sites[1]<<")\n";
		typedef PsimagLite::ContinuedFraction<TridiagonalMatrixType>
		        ContinuedFractionType;
		typedef PsimagLite::ContinuedFractionCollection<ContinuedFractionType>
		        ContinuedFractionCollectionType;

		typename EngineType::VectorStringType vstr;
		PsimagLite::IoSimple::Out ioOut(std::cout);
		ContinuedFractionCollectionType cfCollection(PsimagLite::FREQ_REAL);
		SizeType norbitals = maxOrbitals(model);
		for (SizeType orb1=0;orb1<norbitals;orb1++) {
			for (SizeType orb2=orb1;orb2<norbitals;orb2++) {
				engine.spectralFunction(cfCollection,
				                        vstr,
				                        gfI,
				                        lanczosOptions.sites[0],
				        lanczosOptions.sites[1],
				        lanczosOptions.spins,
				        std::pair<SizeType,SizeType>(orb1,orb2));
			}
		}

		ioOut<<"#INDEXTOCF ";
		for (SizeType i = 0; i < vstr.size(); ++i)
			ioOut<<vstr[i]<<" ";
		ioOut<<"\n";
		cfCollection.save(ioOut);
	}

	for (SizeType cicji=0;cicji<lanczosOptions.cicj.size();cicji++) {
		SizeType cicjI = lanczosOptions.cicj[cicji];
		SizeType total = geometry.numberOfSites();
		PsimagLite::Matrix<ComplexOrRealType> cicjMatrix(total,total);
		SizeType norbitals = maxOrbitals(model);
		for (SizeType orb1=0;orb1<norbitals;orb1++) {
			for (SizeType orb2=0;orb2<norbitals;orb2++) {
				engine.twoPoint(cicjMatrix,
				                cicjI,
				                lanczosOptions.spins,
				                std::pair<SizeType,SizeType>(orb1,orb2));
				std::cout<<cicjMatrix;
			}
		}
	}

	if (lanczosOptions.split >= 0) {
		ReducedDensityMatrix<ModelType> reducedDensityMatrix(model,
		                                                     engine.eigenvector(),
		                                                     lanczosOptions.split);
		reducedDensityMatrix.printAll(std::cout);
	}

	if (lanczosOptions.extendedStatic != "") {
		PsimagLite::Vector<PsimagLite::String>::Type str;
		PsimagLite::tokenizer(lanczosOptions.extendedStatic,str,";");
		VectorSizeType sites;
		VectorSizeType spins;
		VectorSizeType orbs;
		VectorSizeType whats;
		for (SizeType i = 0; i < str.size(); ++i) {
			PsimagLite::Vector<PsimagLite::String>::Type str2;
			PsimagLite::tokenizer(str[i],str2,",");
			if (str2.size() < 3)
				throw PsimagLite::RuntimeError("-S option malformed\n");
			whats.push_back(ProgramGlobals::operator2id(str2[0]));
			sites.push_back(atoi(str2[1].c_str()));
			spins.push_back(atoi(str2[2].c_str()));
			if (str2.size() == 4) orbs.push_back(atoi(str2[3].c_str()));
			else orbs.push_back(0);
		}

		engine.manyPoint(sites,whats,spins,orbs);
	}
}


template<typename ModelType,typename SpecialSymmetryType>
void mainLoop2(const ModelType& model,
               InputNgType::Readable& io,
               LanczosOptions& lanczosOptions)
{
	PsimagLite::String tmp;
	io.readline(tmp,"SolverOptions=");
	bool onthefly = (tmp.find("InternalProductOnTheFly") != PsimagLite::String::npos);

	if (onthefly) {
		mainLoop3<ModelType,SpecialSymmetryType,InternalProductOnTheFly>(model,
		                                                                 io,
		                                                                 lanczosOptions);
	} else {
		mainLoop3<ModelType,SpecialSymmetryType,InternalProductStored>(model,
		                                                               io,
		                                                               lanczosOptions);
	}
}

template<typename ModelType>
void mainLoop(InputNgType::Readable& io,
              const ModelType& model,
              LanczosOptions& lanczosOptions)
{
	typedef typename ModelType::GeometryType GeometryType;
	typedef typename ModelType::BasisBaseType BasisBaseType;

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
		mainLoop2<ModelType,TranslationSymmetry<GeometryType,BasisBaseType> >(model,
		                                                                      io,
		                                                                      lanczosOptions);
	} else if (useReflectionSymmetry) {
		mainLoop2<ModelType,ReflectionSymmetry<GeometryType,BasisBaseType> >(model,
		                                                                     io,
		                                                                     lanczosOptions);
	} else {
		mainLoop2<ModelType,DefaultSymmetry<GeometryType,BasisBaseType> >(model,
		                                                                  io,
		                                                                  lanczosOptions);
	}
}

template<typename ComplexOrRealType>
void mainLoop0(InputNgType::Readable& io,
               LanczosOptions& lanczosOptions)
{
	typedef PsimagLite::Geometry<ComplexOrRealType,
	        InputNgType::Readable,
	        ProgramGlobals> GeometryType;
	typedef ModelSelector<ComplexOrRealType,GeometryType,InputNgType::Readable> ModelSelectorType;
	typedef typename ModelSelectorType::ModelBaseType ModelBaseType;


	GeometryType geometry(io);

	std::cout<<geometry;

	ModelSelectorType modelSelector(io,geometry);
	const ModelBaseType& modelPtr = modelSelector();

	std::cout<<modelPtr;
	mainLoop(io,modelPtr,lanczosOptions);
}

int main(int argc,char *argv[])
{
	int opt = 0;
	LanczosOptions lanczosOptions;
	PsimagLite::String file = "";
	PsimagLite::Vector<PsimagLite::String>::Type str;
	InputCheck inputCheck;
	int precision = 6;
	bool versionOnly = false;

	/* PSIDOC LanczosDriver
	\begin{itemize}
	\item[-g label] Computes the spectral function (continued fraction) for label.
	\item[-c label] Computes the two-point correlation for label.
	\item[-f file] Input file to use. DMRG++ inputs can be used.
	\item[-s ``s1,s2''] computes correlations or spectral functions for spin s1,s2.
	Only s1==s2 is supported for now.
	\item[-r siteForSplit] Calculates the reduced density matrix with a lattice
	split at the siteForSplit.
	\item[-p precision] precision in decimals to use.
	\item[-V] prints version and exits.
	\end{itemize}
	*/
	while ((opt = getopt(argc, argv, "g:c:f:s:r:p:S:V")) != -1) {
		switch (opt) {
		case 'g':
			lanczosOptions.gf.push_back(ProgramGlobals::operator2id(optarg));
			break;
		case 'f':
			file = optarg;
			break;
		case 'c':
			lanczosOptions.cicj.push_back(ProgramGlobals::operator2id(optarg));
			break;
		case 's':
			lanczosOptions.spins.clear();
			PsimagLite::tokenizer(optarg,str,";");
			fillOrbsOrSpin(lanczosOptions.spins,str);
			str.clear();
			break;
		case 'r':
			lanczosOptions.split = atoi(optarg);
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 'S':
			lanczosOptions.extendedStatic = optarg;
			break;
		case 'V':
			versionOnly = true;
			break;
		default: /* '?' */
			inputCheck.usage(argv[0]);
			return 1;
		}
	}

	if (file == "" && !versionOnly) {
		inputCheck.usage(argv[0]);
		return 1;
	}

	//! setup distributed parallelization
	SizeType npthreads = 1;
	ConcurrencyType concurrency(&argc,&argv,npthreads);

	// print license
	if (ConcurrencyType::root()) {
		std::cerr<<license;
		std::cerr<<"Lanczos++ Version "<<LANCZOSPP_VERSION<<"\n";
		std::cerr<<"PsimagLite version "<<PSIMAGLITE_VERSION<<"\n";
	}

	if (versionOnly) return 0;

	//Setup the Geometry
	InputNgType::Writeable ioWriteable(file,inputCheck);
	InputNgType::Readable io(ioWriteable);

	bool isComplex = false;

	PsimagLite::String solverOptions;
	io.readline(solverOptions,"SolverOptions=");
	PsimagLite::Vector<PsimagLite::String>::Type tokens;
	PsimagLite::tokenizer(solverOptions,tokens,",");
	for (SizeType i = 0; i < tokens.size(); ++i) {
		if (tokens[i] == "useComplex") {
			isComplex = true;
			break;
		}
	}

	try {
		io.readline(npthreads,"Threads=");
		ConcurrencyType::npthreads = npthreads;
	} catch (std::exception&) {}

	inputCheck.checkForThreads(ConcurrencyType::npthreads);

	typedef std::complex<RealType> ComplexType;

	if (isComplex)
		mainLoop0<ComplexType>(io,lanczosOptions);
	else
		mainLoop0<RealType>(io,lanczosOptions);
}
