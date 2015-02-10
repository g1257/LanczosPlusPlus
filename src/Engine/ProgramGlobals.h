/*
Copyright (c) 2009-2012, UT-Battelle, LLC
All rights reserved

[LanczosPlusPlus, Version 1.0.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/

/** \ingroup LanczosPlusPlus */
/*@{*/

/*! \file ProgramGlobals.h
 *
 *
 *
 */
#ifndef LANCZOS_PROGRAM_LIMITS_H
#define LANCZOS_PROGRAM_LIMITS_H
#include "String.h"
#include "TypeToString.h"
#include <bitset>
#include <climits>

namespace LanczosPlusPlus {
struct ProgramGlobals {

	typedef std::pair<int,int> PairIntType;
	typedef unsigned int long long WordType;

	enum {FERMION,BOSON};

	enum {SPIN_UP,SPIN_DOWN};

	enum {OPERATOR_NIL,
		  OPERATOR_C,
		  OPERATOR_SZ,
		  OPERATOR_CDAGGER,
		  OPERATOR_N,
		  OPERATOR_SPLUS,
		  OPERATOR_SMINUS};

	static bool needsNewBasis(SizeType what)
	{
		if (what==OPERATOR_C || what==OPERATOR_CDAGGER) return true;
		if (what==OPERATOR_SPLUS || what==OPERATOR_SMINUS) return true;
		return false;
	}

	static SizeType operator2id(const PsimagLite::String& s)
	{
		if (s=="c") {
			return OPERATOR_C;
		} else if (s=="cdagger") {
			return OPERATOR_CDAGGER;
		} else if (s=="sz") {
			return OPERATOR_SZ;
		} else if (s=="nil") {
			return OPERATOR_NIL;
		} else if (s=="n") {
			return OPERATOR_N;
		} else if (s=="splus" || s=="s+") {
			return OPERATOR_SPLUS;
		} else if (s=="sminus" || s=="s-") {
			return OPERATOR_SMINUS;
		}
		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) +  "\n";
		str += "operatorWithType: unsupported operator " + s + "\n";
		throw std::runtime_error(str.c_str());
	}

	static PsimagLite::String id2Operator(SizeType id)
	{
		PsimagLite::Vector<PsimagLite::String>::Type labels;
		labels.push_back("cdagger");
		labels.push_back("c");
		labels.push_back("n");
		labels.push_back("sz");
		labels.push_back("splus");
		labels.push_back("sminus");
		labels.push_back("nil");
		for (SizeType i=0;i<labels.size();i++) {
			if (operator2id(labels[i])==id) return labels[i];
		}
		return "UNKNOWN";
	}

	static PsimagLite::String unknownOperator(SizeType id)
	{
		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Unknown operator " + id2Operator(id) + "\n";
		return str;
	}

	static bool isFermionic(SizeType what)
	{
		if (what==OPERATOR_C) return true;
		if (what==OPERATOR_CDAGGER) return true;
		return false;
	}

	static SizeType transposeConjugate(SizeType operatorLabel)
	{
		if (operatorLabel==OPERATOR_C)
			return OPERATOR_CDAGGER;
		if (operatorLabel==OPERATOR_CDAGGER)
			return OPERATOR_C;
		if (operatorLabel==OPERATOR_SPLUS)
			return OPERATOR_SMINUS;
		if (operatorLabel==OPERATOR_SMINUS)
			return OPERATOR_SPLUS;
		return operatorLabel;
	}

	template<typename T>
	static void binRep(std::ostream& os,
	                   SizeType n,
	                   const T& a)
	{
		const char* beg = reinterpret_cast<const char*>(&a);
		SizeType len = sizeof(a);
		const char* end = beg + len;
		PsimagLite::String buffer("");

		while(beg != end) {
			PsimagLite::String str(std::bitset<CHAR_BIT>(*(end-1)).to_string());
			buffer += str;
			end--;
		}

		SizeType lmn = buffer.length();
		if (lmn >= n) lmn -= n;
		os<<buffer.substr(lmn,n);
	}

	template<typename VectorWordType>
	static void printBasisVector(std::ostream& os,
	                             SizeType n,
	                             VectorWordType& data)
	{
		for (SizeType i=0;i<data.size();i++) {
			binRep(os,n,data[i]);
			os<<"\n";
		}

		os<<"--------------\n";
	}
}; // ProgramGlobals

}; // namespace LanczosPlusPlus
/*@}*/
#endif
