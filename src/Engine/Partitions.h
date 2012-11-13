
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2009-2012, UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. 

Please see full open source license included in file LICENSE.
*********************************************************

*/

#ifndef PARTITIONS_H
#define PARTITIONS_H
#include "Matrix.h"
#include "BitManip.h"

namespace LanczosPlusPlus {

class Partitions {

public:

	Partitions(size_t length,size_t parts)
		: length_(length)
	{
		std::vector<size_t> values(parts,0);

		while(true) {
			if (sumOf(values)==length_)
				partitions_.push_back(values);
			values[0]++;
			if (sumOf(values)>length_) {
				int x = increaseNextIndices(values);
				if (x<0) break;
			}
		}
	}

	size_t size() const { return partitions_.size(); }

	const std::vector<size_t>& operator()(size_t i) const
	{
		return partitions_[i];
	}

private:

	size_t sumOf(const std::vector<size_t>& values) const
	{
		size_t sum = 0;
		for (size_t i=0;i<values.size();i++)
			sum += values[i];
		return sum;
	}

	int increaseNextIndices(std::vector<size_t>& values) const
	{
		if (0==values.size()-1) return -1;
		values[0]=0;
		size_t i = 1;
		while(true) {
			values[i]++;
			if (sumOf(values)<=length_) return i;
			if (i==values.size()-1) return -1;
			values[i]=0;
			i++;
		}
	}

	size_t length_;
	std::vector<std::vector<size_t> > partitions_;

}; // class Partitions



} // namespace LanczosPlusPlus
#endif

