#include "PsimagLite.h"
#include "InputCheck.h"
#include "Concurrency.h"
#include "InputNg.h"
#include "../Version.h"
#include "../../PsimagLite/src/Version.h"

PsimagLite::String license = "Copyright (c) 2009-2021, UT-Battelle, LLC\n"
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

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif
typedef std::complex<RealType> ComplexType;
typedef PsimagLite::InputNg<InputCheck> InputNgType;

template<typename ComplexOrRealType>
void mainLoop0(InputNgType::Readable& io)
{

}

int main(int argc,char **argv)
{
	PsimagLite::PsiApp application("lanczos++",&argc,&argv,1);
	PsimagLite::String file;
	int precision = 0;
	bool versionOnly = false;
	SizeType threadsInCmdLine = 0;
	InputCheck inputCheck;

	/* PSIDOC Dynamics1Driver
	\begin{itemize}
	\item[-f file] Input file to use. DMRG++ inputs can be used.
	\item[-S number] Override Threads= in input line if preset, and set threads
	to number given here.
	\item[-p precision] precision in decimals to use.
	\item[-V] prints version and exits.
	\end{itemize}
	*/
	int opt = 0;
	while ((opt = getopt(argc, argv, "g:c:m:f:s:r:p:M:S:V")) != -1) {
		switch (opt) {
		case 'f':
			file = optarg;
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 'S':
			threadsInCmdLine = atoi(optarg);
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

	// print license
	if (PsimagLite::Concurrency::root()) {
		std::cerr<<license;
		std::cerr<<"Lanczos++ Version "<<LANCZOSPP_VERSION<<"\n";
		std::cerr<<"PsimagLite version "<<PSIMAGLITE_VERSION<<"\n";
	}

	if (versionOnly) return 0;

	InputNgType::Writeable ioWriteable(file,inputCheck);
	InputNgType::Readable io(ioWriteable);

	bool isComplex = false;
	bool setAffinities = false;

	PsimagLite::String solverOptions;
	io.readline(solverOptions,"SolverOptions=");
	PsimagLite::Vector<PsimagLite::String>::Type tokens;
	PsimagLite::split(tokens, solverOptions, ",");
	for (SizeType i = 0; i < tokens.size(); ++i) {
		if (tokens[i] == "useComplex") {
			isComplex = true;
		} else if (tokens[i] == "setAffinities") {
			setAffinities = true;
		}
	}

	//! setup distributed parallelization
	SizeType npthreads = 1;
	try {
		io.readline(npthreads,"Threads=");
	} catch (std::exception&) {}

	if (threadsInCmdLine > 0)
		npthreads = threadsInCmdLine;

	PsimagLite::CodeSectionParams codeSectionParams(npthreads, 1, setAffinities, 0);
	PsimagLite::Concurrency::setOptions(codeSectionParams);

	typedef std::complex<RealType> ComplexType;

	if (isComplex)
		mainLoop0<ComplexType>(io);
	else
		mainLoop0<RealType>(io);
}

