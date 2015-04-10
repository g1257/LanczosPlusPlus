#include <iostream>
#include <fstream>
#include <cassert>
#include <unistd.h>
#include "Vector.h"
#include "Sort.h"

typedef double RealType;
typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

void load(VectorRealType& e, VectorRealType& w,PsimagLite::String file)
{
	std::ifstream fin(file.c_str());
	RealType tmp = 0;
	while (!fin.eof() && !fin.bad() && fin.good()) {
		fin>>tmp;
		e.push_back(tmp);
		fin>>tmp;
		w.push_back(tmp);
	}

	std::cerr<<"load: "<<e.size()<<" values found in "<<file<<"\n";
}

void sort(VectorRealType& e, VectorRealType& w)
{
	assert(e.size() == w.size());
	if (e.size() == 0) return;

	PsimagLite::Sort<VectorRealType> sort;
	VectorSizeType iperm(e.size());
	sort.sort(e,iperm);
	VectorRealType ww(w.size());
	for (SizeType i = 0; i < w.size(); ++i) {
		ww[i] = w[iperm[i]];
	}

	w = ww;
}

void prune(VectorRealType& e,
           VectorRealType& w,
           RealType& emin,
           RealType& emax,
           RealType& wabsmax)
{
	assert(e.size() == w.size());
	if (e.size() == 0) return;

	SizeType i = 0;
	for (; i < w.size(); ++i) {
		if (fabs(w[i])>1e-6) break;
	}

	if (i > 0) i--;

	int j = w.size() - 1;
	for (; j >= 0; j--) {
		if (fabs(w[j])>1e-6) break;
	}

	SizeType final = static_cast<SizeType>(j+1);

	SizeType counter = 0;
	SizeType total = (final - i);
	VectorRealType ee(total);
	VectorRealType ww(total);
	for (SizeType k = i; k < final; ++k) {
		ee[counter] = e[k];
		ww[counter++] = w[k];
		if (e[k] > emax) emax = e[k];
		if (e[k] < emin) emin = e[k];
		if (fabs(w[k]) > wabsmax) wabsmax = fabs(w[k]);
	}

	e = ee;
	w = ww;
	std::cerr<<"prun: emax="<<emax<<" emin="<<emin<<" wabsmax="<<wabsmax<<"\n";
	std::cerr<<"prune: "<<e.size()<<" values remain after pruning\n";
}

RealType lorentzian(RealType omega,
                    const VectorRealType& e,
                    const VectorRealType& w,
                    RealType eps)
{
	RealType eps2 = eps*eps;
	RealType sum = 0;
	for (SizeType i = 0; i < e.size(); ++i) {
		RealType tmp = (omega - e[i]);
		sum += w[i]/(tmp*tmp+eps2);
	}

	return sum;
}

void usage(char *name, PsimagLite::String msg = "")
{
	if (msg != "") std::cerr<<name<<": "<<msg<<"\n";
	std::cerr<<"USAGE: "<<name<<" -f file -e eps -t total\n";
}

int main(int argc, char **argv)
{
	PsimagLite::String file;
	RealType eps = 0.1;
	SizeType total = 0;
	int opt = 0;
	while ((opt = getopt(argc, argv, "f:e:t:")) != -1) {
		switch (opt) {
		case 'f':
			file = optarg;
			break;
		case 'e':
			eps = atof(optarg);
			break;
		case 't':
			total = atoi(optarg);
			break;
		default: /* '?' */
			usage(argv[0]);
			return 1;
		}
	}

	if (file == "" || total == 0) {
		usage(argv[0]);
		return 2;
	}

	// load
	VectorRealType e;
	VectorRealType w;
	load(e,w,file);
	// sort
	sort(e,w);
	// min, max, prune
	RealType emin = 1e10;
	RealType emax = -emin;
	RealType wabsmax = 0;
	prune(e,w,emin,emax,wabsmax);

	RealType omegaInit = emin; // FIXME: allow override
	RealType omegaStep = (emax-omegaInit)/(total-1); // FIXME: allow override
	RealType factor = eps/wabsmax;

	for (SizeType i = 0; i < total; ++i) {
		RealType omega = i*omegaStep + omegaInit;
		RealType val = lorentzian(omega,e,w,eps);
		val *= factor;
		std::cout<<omega<<" "<<val<<"\n";
	}
}

