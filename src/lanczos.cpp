#include "LanczosDriver.h"

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

using namespace LanczosPlusPlus;

void fillOrbsOrSpin(PsimagLite::Vector<LanczosPlusPlus::LanczosOptions::PairSizeType>::Type& spinV,
                    const PsimagLite::Vector<PsimagLite::String>::Type& strV)
{
	for (SizeType i=0;i<strV.size();i++) {
		PsimagLite::Vector<PsimagLite::String>::Type strV2;
		PsimagLite::tokenizer(strV[i],strV2,",");
		if (strV2.size()!=2)
			throw std::runtime_error("-o needs pairs\n");
		LanczosPlusPlus::LanczosOptions::PairSizeType spins;
		spins.first = atoi(strV2[0].c_str());
		spins.second = atoi(strV2[1].c_str());
		spinV.push_back(spins);
	}
}

template<typename ModelType>
void mainLoop(InputNgType::Readable& io,
              const ModelType& model,
              LanczosPlusPlus::LanczosOptions& lanczosOptions)
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
		mainLoop2<ModelType,
		        LanczosPlusPlus::TranslationSymmetry<GeometryType,BasisBaseType> >(model,
		                                                                           io,
		                                                                           lanczosOptions);
	} else if (useReflectionSymmetry) {
		mainLoop2<ModelType,
		        LanczosPlusPlus::ReflectionSymmetry<GeometryType,BasisBaseType> >(model,
		                                                                          io,
		                                                                          lanczosOptions);
	} else {
		mainLoop2<ModelType,
		        LanczosPlusPlus::DefaultSymmetry<GeometryType,BasisBaseType> >(model,
		                                                                       io,
		                                                                       lanczosOptions);
	}
}

template<typename ComplexOrRealType>
void mainLoop0(InputNgType::Readable& io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions)
{
	typedef PsimagLite::Geometry<ComplexOrRealType,
	        InputNgType::Readable,
	        LanczosPlusPlus::ProgramGlobals> GeometryType;
	typedef LanczosPlusPlus::ModelSelector<ComplexOrRealType,
	        GeometryType,
	        InputNgType::Readable> ModelSelectorType;
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
	PsimagLite::PsiApp application("lanczos++");
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
