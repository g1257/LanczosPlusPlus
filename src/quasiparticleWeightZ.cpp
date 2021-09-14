#include "lanczos.cpp"
#include "LabeledOperator.h"


template<typename EngineType, bool b>
class FindVectorTtoGs {

public:

	typedef typename EngineType::ComplexOrRealType ComplexOrRealType;
	typedef typename EngineType::ModelType ModelType;

	void findVectorTtoGs(std::vector<ComplexOrRealType>& zvector,
	                     const EngineType& engine,
	                     const ModelType& model,
	                     SizeType mForK)
	{
		throw PsimagLite::RuntimeError("Only works with useComplex\n");
	}
};

template<typename EngineType>
class FindVectorTtoGs<EngineType, true> {

public:

	typedef typename EngineType::ComplexOrRealType ComplexOrRealType;
	typedef typename EngineType::ModelType ModelType;
	typedef std::pair<SizeType, SizeType> PairSizeType;

	FindVectorTtoGs() : basisNew_(nullptr) {}

	void findVectorTtoGs(std::vector<ComplexOrRealType>& zvector,
	                     const EngineType& engine,
	                     const ModelType& model,
	                     SizeType mForK)
	{
		const SizeType sites = model.geometry().numberOfSites();
		const std::vector<ComplexOrRealType>& gsVector = engine.eigenvector(0);
		const SizeType spin = 0; // bogus
		const SizeType orb = 0; // bogus
		const PairSizeType spins(spin, spin);
		static const LabeledOperator::Label label = LabeledOperator::Label::OPERATOR_CDAGGER;

		createBasisIfNeeded(label, model, spins, orb);

		const SizeType n = (basisNew_) ? basisNew_->size() : gsVector.size();
		zvector.resize(n);

		for (SizeType site = 0; site < sites; ++site) {
			RealType arg = 2.0*M_PI*mForK/sites;
			ComplexOrRealType factor(cos(arg), sin(arg));
			engine.accModifiedState_(zvector,
			                         label,
			                         (basisNew_) ? *basisNew_ : model.basis(),
			                         gsVector,
			                         model.basis(),
			                         site,
			                         spin,
			                         orb,
			                         factor);
		}
	}

private:

	void createBasisIfNeeded(const LabeledOperator& lOperator,
	                         const ModelType& model,
	                         const PairSizeType& spins,
	                         SizeType orb)
	{
		if (lOperator.needsNewBasis()) {
			if (spins.first != spins.second) {
				PsimagLite::String str(__FILE__);
				str += " " + ttos(__LINE__) + "\n";
				str += "twoPoint: no support yet for off-diagonal spin ";
				str += "when needs new basis\n";
				throw std::runtime_error(str.c_str());
			}

			PairSizeType oldParts = model.basis().parts();
			PairSizeType newParts(0,0);
			if (!model.hasNewParts(newParts, oldParts, lOperator, spins.first, orb))
				return;

			basisNew_ = model.createBasis(newParts.first, newParts.second);

			std::cerr<<"basisNew.size="<<basisNew_->size()<<" ";
			std::cerr<<"newparts.first="<<newParts.first<<" ";
			std::cerr<<"newparts.second="<<newParts.second<<"\n";
		}
	}

	typename EngineType::BasisType* basisNew_;
};

bool discardThisLine(PsimagLite::String line)
{
    PsimagLite::Vector<PsimagLite::String>::Type lines = {"TargetElectronsUp=",
                                                          "TargetElectronsDown="};
    const SizeType n = lines.size();
    for (SizeType i = 0; i < n; ++i) {
        std::size_t found = line.find(lines[i]);
        if (found != PsimagLite::String::npos)
            return true;
    }

    return false;
}

PsimagLite::String makeDataForOneHole(const PsimagLite::String& str,
                                      SizeType nup,
                                      SizeType ndown,
                                      bool isAinur)
{
    const SizeType total = str.length();
    PsimagLite::String buffer;
    PsimagLite::String str2;
    for (SizeType i = 0; i < total; ++i) {
        const char c = str[i];
        buffer.push_back(c);
        if (c == '\n') {
            if (!discardThisLine(buffer))
                str2 += buffer;
            buffer = "";
        }
    }

    PsimagLite::String semiOrNothing = (isAinur) ? ";" : "";
    str2 += "TargetElectronsUp=" + ttos(nup) + semiOrNothing + "\n";
    str2 += "TargetElectronsDown=" + ttos(ndown) + semiOrNothing + "\n";;
    return str2;
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
	typedef LanczosPlusPlus::ModelSelector<ComplexOrRealType,
	        GeometryType,
	        InputNgType::Readable> ModelSelectorType;
	typedef std::pair<SizeType, SizeType> PairSizeType;

	EngineType engine(model, io);

	//! get the g.s.:
	RealType Eg = engine.energies(0);
	std::cout.precision(8);
	std::cout<<"Energy="<<Eg<<"\n";

	std::vector<ComplexOrRealType> tTogs;
	FindVectorTtoGs<EngineType, PsimagLite::IsComplexNumber<ComplexOrRealType>::True> fv;
	fv.findVectorTtoGs(tTogs, engine, model, lanczosOptions.split);

	ComplexOrRealType den = scalarProduct(tTogs, tTogs);
	std::cout<<"Denominator="<<den<<"\n";

	PairSizeType oldParts = model.basis().parts();

	InputCheck inputCheck;
	PsimagLite::String input0 = makeDataForOneHole(io.data(),
	                                               oldParts.first - 1, // nup - 1
	                                               oldParts.second, // ndown
	                                               io.isAinur());
	//std::cout<<input0;
	InputNgType::Writeable ioWriteable2(inputCheck, input0);
	InputNgType::Readable ioOneHole(ioWriteable2);

	ModelSelectorType modelSelector2(ioOneHole, model.geometry());

	EngineType engine2(modelSelector2(), ioOneHole);
	RealType EgOneHole = engine2.energies(0);
	std::cout.precision(8);
	std::cout<<"EnergyOneHole="<<EgOneHole<<"\n";
	ComplexOrRealType num = scalarProduct(engine2.eigenvector(0), tTogs);
	std::cout<<"Numerator="<<num<<"\n";
	std::cout<<"QuasiparticleWeightZ="<<num/den<<"\n";
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
