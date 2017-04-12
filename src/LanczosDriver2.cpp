#include "LanczosDriver1.h"

typedef PsimagLite::Geometry<RealType,
        InputNgType::Readable,
        LanczosPlusPlus::ProgramGlobals> Geometry2Type;

typedef LanczosPlusPlus::ModelSelector<RealType,
        Geometry2Type,
        InputNgType::Readable> ModelSelector2Type;
typedef typename ModelSelector2Type::ModelBaseType ModelBase2Type;

typedef typename ModelBase2Type::BasisBaseType BasisBase2Type;

typedef LanczosPlusPlus::DefaultSymmetry<Geometry2Type,BasisBase2Type> Symmetry2Type;

template
void mainLoop2<ModelBase2Type,Symmetry2Type>(const ModelBase2Type& model,
               InputNgType::Readable& io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions);
