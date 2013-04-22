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
#ifndef PROGRAM_LIMITS_H
#define PROGRAM_LIMITS_H
#include <string>
#include "TypeToString.h"

namespace LanczosPlusPlus {
struct ProgramGlobals {
	static size_t const MaxLanczosSteps = 1000000; // max number of internal Lanczos steps
	static size_t const LanczosSteps = 300; // max number of external Lanczos steps
	static double const LanczosTolerance; // tolerance of the Lanczos Algorithm

	enum {FERMION,BOSON};

	enum {OPERATOR_NIL,
	      OPERATOR_C,
	      OPERATOR_SZ,
	      OPERATOR_CDAGGER,
	      OPERATOR_N,
	      OPERATOR_SPLUS,
	      OPERATOR_SMINUS};

	static bool needsNewBasis(size_t what)
	{
		if (what==OPERATOR_C) return true;
		if (what==OPERATOR_CDAGGER) return true;
		return false;
	}

	static size_t operator2id(const std::string& s)
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
		} else if (s=="splus") {
			return OPERATOR_SPLUS;
		} else if (s=="sminus") {
			return OPERATOR_SMINUS;
		}
		std::string str = unknownOperator(s);
		throw std::runtime_error(str.c_str());
	}

	static std::string id2Operator(size_t id)
	{
		std::vector<std::string> labels;
		labels.push_back("cdagger");
		labels.push_back("c");
		labels.push_back("n");
		labels.push_back("sz");
		labels.push_back("splus");
		labels.push_back("sminus");
		labels.push_back("nil");
		for (size_t i=0;i<labels.size();i++) {
			if (operator2id(labels[i])==id) return labels[i];
		}
		return "UNKNOWN";
	}

	static size_t operatorWithType(size_t op,size_t type)
	{
		if (op==OPERATOR_C || op==OPERATOR_CDAGGER) {
			if (type&1) return transposeConjugate(op);
			return op;
		}
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) +  "\n";
		str += "operatorWithType: unsupported operator " + id2Operator(op) + "\n";
		throw std::runtime_error(str.c_str());

	}

	static std::string unknownOperator(const std::string& s)
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Unknown operator " + s + "\n";
		return str;
	}

	static bool isFermionic(size_t what)
	{
		if (what==OPERATOR_C) return true;
		if (what==OPERATOR_CDAGGER) return true;
		return false;
	}

	static size_t transposeConjugate(size_t operatorLabel)
	{
		if (operatorLabel==OPERATOR_C)
			return OPERATOR_CDAGGER;
		if (operatorLabel==OPERATOR_CDAGGER)
			return OPERATOR_C;
		return operatorLabel;
	}
}; // ProgramGlobals

	double const ProgramGlobals::LanczosTolerance = 1e-12;
}; // namespace LanczosPlusPlus
/*@}*/
#endif
