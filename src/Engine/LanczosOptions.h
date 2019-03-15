#ifndef LANCZOSOPTIONS_H
#define LANCZOSOPTIONS_H
#include "Vector.h"
#include "LabeledOperator.h"

namespace  LanczosPlusPlus {

struct LanczosOptions {

	typedef std::pair<SizeType,SizeType> PairSizeType;

	LanczosOptions()
	    : split(-1),spins(1,PairSizeType(0,0)),extendedStatic("")
	{}

	int split;
	PsimagLite::Vector<LabeledOperator>::Type cicj;
	PsimagLite::Vector<LabeledOperator>::Type gf;
	PsimagLite::Vector<SizeType>::Type sites;
	PsimagLite::Vector<PairSizeType>::Type spins;
	PsimagLite::String extendedStatic;

}; // struct LanczosOptions
}

#endif // LANCZOSOPTIONS_H
