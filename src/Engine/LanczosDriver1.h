#ifndef LANCZOSDRIVER_1_H
#define LANCZOSDRIVER_1_H
#include "LanczosDriver.h"

template<typename ModelType>
SizeType maxOrbitals(const ModelType& model)
{
	SizeType res=0;
	for (SizeType i=0;i<model.geometry().numberOfSites();i++) {
		if (res<model.orbitals(i)) res=model.orbitals(i);
	}
	return res;
}

template<typename EngineType>
void extendedStatic(PsimagLite::String manypoint, const EngineType& engine)
{
	typedef typename EngineType::VectorSizeType VectorSizeType;
	PsimagLite::Vector<PsimagLite::String>::Type str;
	PsimagLite::split(str, manypoint, ";");

	VectorSizeType sites;
	VectorSizeType spins;
	VectorSizeType orbs;
	VectorSizeType whats;
	for (SizeType i = 0; i < str.size(); ++i) {
		PsimagLite::Vector<PsimagLite::String>::Type str2;
		PsimagLite::split(str2, str[i], "?");
		if (str2.size() < 3)
			throw PsimagLite::RuntimeError("-S option malformed\n");
		whats.push_back(LanczosPlusPlus::ProgramGlobals::operator2id(str2[0]));
		sites.push_back(atoi(str2[1].c_str()));
		spins.push_back(atoi(str2[2].c_str()));
		if (str2.size() == 4) orbs.push_back(atoi(str2[3].c_str()));
		else orbs.push_back(0);
	}

	std::cout<<"<gs|"<<manypoint<<"|gs>=";
	std::cout<<engine.manyPoint(sites,whats,spins,orbs);
	std::cout<<"\n";
}

template<typename ModelType,
         typename SpecialSymmetryType,
         template<typename,typename> class InternalProductTemplate>
void mainLoop3(const ModelType& model,
               InputNgType::Readable& io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions)
{
	typedef typename ModelType::GeometryType GeometryType;
	typedef typename GeometryType::ComplexOrRealType ComplexOrRealType;
	typedef LanczosPlusPlus::Engine<ModelType,InternalProductTemplate,SpecialSymmetryType> EngineType;
	typedef typename EngineType::TridiagonalMatrixType TridiagonalMatrixType;

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
		LanczosPlusPlus::ReducedDensityMatrix<ModelType> reducedDensityMatrix(model,
		                                                                      engine.eigenvector(),
		                                                                      lanczosOptions.split);
		reducedDensityMatrix.printAll(std::cout);
	}

	if (lanczosOptions.extendedStatic != "") {
		PsimagLite::Vector<PsimagLite::String>::Type str;
		PsimagLite::split(str, lanczosOptions.extendedStatic, ",");
		for (SizeType i = 0; i < str.size(); ++i)
			extendedStatic(str[i],engine);
	}
}


template<typename ModelType,typename SpecialSymmetryType>
void mainLoop2(const ModelType& model,
               InputNgType::Readable& io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions)
{
	PsimagLite::String tmp;
	io.readline(tmp,"SolverOptions=");
	bool onthefly = (tmp.find("InternalProductOnTheFly") != PsimagLite::String::npos);

	if (onthefly) {
		mainLoop3<ModelType,
		        SpecialSymmetryType,
		        LanczosPlusPlus::InternalProductOnTheFly>(model,
		                                                  io,
		                                                  lanczosOptions);
	} else {
		mainLoop3<ModelType,
		        SpecialSymmetryType,
		        LanczosPlusPlus::InternalProductStored>(model,
		                                                io,
		                                                lanczosOptions);
	}
}

#endif // LANCZOSDRIVER_1_H
