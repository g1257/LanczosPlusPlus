#ifndef LANCZOSDRIVER_H
#define LANCZOSDRIVER_H
#include "AllocatorCpu.h"
#include "PsimagLite.h"
#include "../Version.h"
#include "../../../PsimagLite/src/Version.h"
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
#include "InputCheck.h"
#include "ReducedDensityMatrix.h"
#include "LanczosOptions.h"

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif

typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::InputNg<LanczosPlusPlus::InputCheck> InputNgType;

template<typename ModelType>
SizeType maxOrbitals(const ModelType& model);

template<typename EngineType>
void extendedStatic(PsimagLite::String manypoint,
                    const EngineType& engine,
                    const typename EngineType::PairType& braAndKet);

template<typename ModelType,
         typename SpecialSymmetryType,
         template<typename,typename> class InternalProductTemplate>
void mainLoop3(const ModelType& model,
               InputNgType::Readable& io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions);


template<typename ModelType,typename SpecialSymmetryType>
void mainLoop2(const ModelType& model,
               InputNgType::Readable& io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions);

template<typename ModelType>
void mainLoop(InputNgType::Readable& io,
              const ModelType& model,
              LanczosPlusPlus::LanczosOptions& lanczosOptions);

template<typename ComplexOrRealType>
void mainLoop0(InputNgType::Readable& io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions);

#endif // LANCZOSDRIVER_H
