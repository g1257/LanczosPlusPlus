#include "lanczos.cpp"
#include "LabeledOperator.h"

template<typename EngineType, typename ComplexOrRealType, typename ModelType>
void findVectorTtoGs(std::vector<ComplexOrRealType>& zvector,
                     const EngineType& engine,
                     const ModelType& model)
{
	const SizeType sites = model.geometry().numberOfSites();
	const std::vector<ComplexOrRealType>& gsVector = engine.eigenvector(0);
	const SizeType n = gsVector.size();
	const SizeType spin = 0; // bogus
	const SizeType orb = 0; // bogus
	const RealType isign = 1;
	zvector.resize(n);
	for (SizeType site = 0; site < sites; ++site) {
		engine.accModifiedState_(zvector,
		                         LabeledOperator::Label::OPERATOR_CDAGGER_A_UP_C_B_UP,
		                         model.basis(),
		                         gsVector,
		                         model.basis(),
		                         site,
		                         spin,
		                         orb,
		                         isign);
	}
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
	//typedef typename EngineType::TridiagonalMatrixType TridiagonalMatrixType;
	//typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	EngineType engine(model, io);

	//! get the g.s.:
	RealType Eg = engine.energies(0);
	std::cout.precision(8);
	std::cout<<"Energy="<<Eg<<"\n";

	std::vector<ComplexOrRealType> tTogs;
	findVectorTtoGs(tTogs, engine, model);
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
