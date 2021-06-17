#include "lanczos.cpp"


template<typename ModelType,
         typename SpecialSymmetryType,
         template<typename,typename> class InternalProductTemplate>
void mainLoop3(const ModelType& model,
               InputNgType::Readable& io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions)
{
	//typedef typename ModelType::GeometryType GeometryType;
	//typedef typename GeometryType::ComplexOrRealType ComplexOrRealType;
	typedef LanczosPlusPlus::Engine<ModelType,InternalProductTemplate,SpecialSymmetryType> EngineType;
	//typedef typename EngineType::TridiagonalMatrixType TridiagonalMatrixType;
	//typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	//const GeometryType& geometry = model.geometry();
	EngineType engine(model, io);

	//! get the g.s.:
	RealType Eg = engine.energies(0);
	std::cout.precision(8);
	std::cout<<"Energy="<<Eg<<"\n";
	std::cerr<<"NOT DONE YET\n";
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
