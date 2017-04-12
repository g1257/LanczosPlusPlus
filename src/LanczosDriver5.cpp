#include "LanczosDriver1.h"

typedef PsimagLite::Geometry<std::complex<RealType>,
        InputNgType::Readable,
        LanczosPlusPlus::ProgramGlobals> Geometry5Type;

typedef LanczosPlusPlus::ModelSelector<std::complex<RealType>,
        Geometry5Type,
        InputNgType::Readable> ModelSelector5Type;
typedef typename ModelSelector5Type::ModelBaseType ModelBase5Type;

typedef typename ModelBase5Type::BasisBaseType BasisBase5Type;

typedef LanczosPlusPlus::DefaultSymmetry<Geometry5Type,BasisBase5Type> Symmetry5Type;

template
void mainLoop2<ModelBase5Type,Symmetry5Type>(const ModelBase5Type& model,
               InputNgType::Readable& io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions);
