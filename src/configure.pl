#!/usr/bin/perl 
=pod
// BEGIN LICENSE BLOCK
Copyright (c) 2009 , UT-Battelle, LLC
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


// END LICENSE BLOCK
=cut
use warnings;
use strict;

my $hasGsl = "no"; # say "no" here to remove GSL dependence

my $mpi=0;
my $platform="linux";
my $lapack="-llapack";
my $PsimagLite="../../PsimagLite/src";
my ($pthreads,$pthreadsLib)=(0,"");
my $brand= "v1.0";


my $gslLibs = " -lgsl  -lgslcblas ";
$gslLibs =" " if ($hasGsl=~/n/i);

#guessPlatform();

#welcome();

my $model="";
my $modelLocation="";
my $stored="";
askQuestions();

createMakefile();

createDriver();


sub askQuestions
{
	print "Enter model\n";
	print "Available: HubbardOneOrbital FeBasedSc\n";
	print "Default: HubbardOneOrbital (press ENTER): ";
	$_=<STDIN>;
	chomp;
	s/ //g;
	if ($_ eq "") {
		$_="HubbardOneOrbital";
	}
	$model = $_;
	$modelLocation = "-IModels/$model";

	print "Do you want matrix-vector to be Stored or OnTheFly?\n";
	print "Available: Stored OnTheFly\n";
	print "Default: Stored (press ENTER): ";
	$_=<STDIN>;
	chomp;
	s/ //g;
	if ($_ eq "") {
		$_="Stored";
	}
	$stored = $_;
}

sub createMakefile
{
	system("cp Makefile Makefile.bak") if (-r "Makefile");
	my $compiler = compilerName();
	open(FOUT,">Makefile") or die "Cannot open Makefile for writing: $!\n";
	my $usePthreadsOrNot = " ";
	$usePthreadsOrNot = " -DUSE_PTHREADS " if ($pthreads);
	if ($modelLocation=~/extendedhubbard1orb/i) {
		$modelLocation = $modelLocation." -IModels/HubbardOneBand ";
	}

print FOUT<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl
# Lanczos++ ($brand) by G.A.
# Platform: $platform
# MPI: $mpi

LDFLAGS =    $lapack  $gslLibs $pthreadsLib
CPPFLAGS = -Werror -Wall  -IEngine $modelLocation -I$PsimagLite/Geometry -I$PsimagLite $usePthreadsOrNot
EOF
if ($mpi) {
	print FOUT "CXX = mpicxx -O3 -DNDEBUG \n";
} else {
	print FOUT "CXX = $compiler  -O3 -DNDEBUG\n";
	print FOUT "#Comment out line below for debugging: \n";
	print FOUT "#CXX = $compiler -g3 \n";
	print FOUT "#Comment out line below for valgrind callgrind analysis: \n";
	print FOUT "#CXX = $compiler -g3 -DNDEBUG -O1 \n";
}
print FOUT<<EOF;
EXENAME = lanczos
all: \$(EXENAME)

lanczos.cpp: configure.pl
	perl configure.pl

lanczos:  lanczos.o 
	\$(CXX) -o lanczos lanczos.o \$(LDFLAGS)  

# dependencies brought about by Makefile.dep
%.o: %.cpp Makefile
	\$(CXX) \$(CPPFLAGS) -c \$< 

Makefile.dep: lanczos.cpp
	\$(CXX) \$(CPPFLAGS) -MM lanczos.cpp  > Makefile.dep

clean:
	rm -f core* \$(EXENAME) *.o

include Makefile.dep

######## End of Makefile ########

EOF
	close(FOUT);
	print STDERR "File Makefile has been written\n";
}


sub computeBackwardMovements
{
	my ($directory) = @_;
	my $ret = "";
	while ($directory =~ s/\///) {
		$ret = $ret."../";
	}
	return $ret;	
}

sub createDriver
{
	my $driverName = "lanczos";
	system("cp $driverName.cpp $driverName.bak") if (-r "$driverName.cpp");
	open(FOUT,">$driverName.cpp") or die "Cannot open file $driverName.cpp for writing: $!\n";
	my $license=getLicense();
	my $concurrencyName = "ConcurrencySerial"; #getConcurrencyName();
	my $parametersName = getParametersName();
	my $pthreadsName = "NoPthreads.h"; #getPthreadsName();
	my $operatorsName = getOperatorsName();
	
print FOUT<<EOF;
/* DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
 * This driver program was written by configure.pl
 * Lanczos++ ($brand) by G.A.*/
#include <unistd.h>
#include <cstdlib>
#include <getopt.h>
#include "$concurrencyName.h"
#include "Engine.h"
#include "$model.h"
#include "$parametersName.h"
#include "Geometry.h"
#include "InternalProduct$stored.h"
#include "IoSimple.h" // in PsimagLite
#include "ProgramGlobals.h"
#include "ContinuedFraction.h" // in PsimagLite 
#include "ContinuedFractionCollection.h" // in PsimagLite
#include "ReflectionSymmetry.h"
#include "Split.h"

using namespace LanczosPlusPlus;

typedef double RealType;
typedef std::complex<RealType> ComplexType;
typedef PsimagLite::$concurrencyName<RealType> ConcurrencyType;
typedef PsimagLite::Geometry<RealType,ProgramGlobals> GeometryType;
typedef $parametersName<RealType> ParametersModelType;
typedef PsimagLite::IoSimple::In IoInputType;
typedef $model<RealType,ParametersModelType,GeometryType> ModelType;
typedef ModelType::BasisType BasisType;
typedef ReflectionSymmetry<GeometryType,BasisType> ReflectionSymmetryType;
typedef InternalProduct$stored<ModelType,ReflectionSymmetryType> InternalProductType;
typedef Engine<ModelType,InternalProductType,ConcurrencyType> EngineType;
typedef EngineType::TridiagonalMatrixType TridiagonalMatrixType;

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" [-g -c] -f filename\\n";
}

int main(int argc,char *argv[])
{
	int opt = 0;
	bool gf = false;
	std::string file = "";
	std::vector<size_t> sites;
	bool cicj=false;
	while ((opt = getopt(argc, argv, "gcf:")) != -1) {
		switch (opt) {
		case 'g':
			gf = true;
			break;
		case 'f':
			file = optarg;
			break;
		case 'c':
			cicj = true;
			break;
		default: /* '?' */
			usage(argv[0]);
			return 1;
		}
	}
	if (file == "") {
		usage(argv[0]);
		return 1;
	}
	//! setup distributed parallelization
	ConcurrencyType concurrency(argc,argv);

	//Setup the Geometry
	IoInputType io(file);
	GeometryType geometry(io);

	// read model parameters
	ParametersModelType mp(io);

	// print license
	std::string license = $license;
	if (concurrency.root()) std::cerr<<license;

	std::vector<RealType> qns;
	io.read(qns,"TargetQuantumNumbers");
	if (qns.size()<2) throw std::runtime_error("HubbardLanczos::ctor(...)\\n");
	size_t nup=size_t(geometry.numberOfSites()*qns[0]);
	size_t ndown=size_t(geometry.numberOfSites()*qns[1]);

	//! Setup the Model
	ModelType model(nup,ndown,mp,geometry);

	EngineType engine(model,geometry.numberOfSites(),io);

	//! get the g.s.:
	RealType Eg = engine.gsEnergy();
	std::cout<<"Energy="<<Eg<<"\\n";
	if (gf) {
		io.read(sites,"TSPSites");
		if (sites.size()==0) throw std::runtime_error("No sites in input file!\\n");
		if (sites.size()==1) sites.push_back(sites[0]);

		std::cout<<"#gf(i="<<sites[0]<<",j="<<sites[1]<<")\\n";
		typedef PsimagLite::ContinuedFraction<RealType,TridiagonalMatrixType>
		ContinuedFractionType;
		typedef PsimagLite::ContinuedFractionCollection<ContinuedFractionType>
			ContinuedFractionCollectionType;

		ContinuedFractionCollectionType cfCollection;
		engine.greenFunction(cfCollection,sites[0],sites[1],ModelType::SPIN_UP,std::pair<size_t,size_t>(0,0));
	
	
		PsimagLite::IoSimple::Out ioOut(std::cout);
		cfCollection.save(ioOut);
	}
	if (!cicj) return 0;

	size_t total = geometry.numberOfSites();
	PsimagLite::Matrix<RealType> cicjMatrix(total,total);
	size_t norbitals = model.orbitals();
	for (size_t i=0;i<norbitals;i++) {
		engine.ciCj(cicjMatrix,ModelType::SPIN_UP,std::pair<size_t,size_t>(i,i));
		std::cout<<cicjMatrix;
	}
}

EOF
	close(FOUT);
	print STDERR "File dmrg.cpp has been written\n";

}



sub getLicense
{
	open(THISFILE,"$0") or return " ";
	while(<THISFILE>) {
		last if (/BEGIN LICENSE BLOCK/);
	}
	my $l="";
	while(<THISFILE>) {
		chomp;
		s/\"/\\\"/g;
		$l = "$l\"$_\\n\"\n";
		last if (/END LICENSE BLOCK/);
	}
	close(THISFILE);
	return $l;
}


sub guessPlatform
{
	$platform="Linux";
	$platform="Darwin" if (isAMac());
}

sub isAMac
{
	open(PIPE,"uname -a |grep -i Darwin") or return 0;
	$_=<PIPE>;
	close(PIPE);
	
	return 1 unless ($_ eq "" or $_ eq "\n");
	return 0;
}


sub getParametersName
{
	my $parametersName="Parameters$model";
	if ($model=~/hubbard/i) {
		$parametersName = "ParametersModelHubbard";
	} elsif ($model=~/heisenberg/i) {
		$parametersName = "ParametersModelHeisenberg";
	} elsif ($model=~/febasedsc/i) {
		$parametersName = "ParametersModelFeAs";
	} elsif ($model=~/tjoneorbital/i) {
		$parametersName = "ParametersTjOneOrbital";
	}
	return $parametersName;
}

sub getPthreadsName
{
	my $pthreadsName = "UNKNOWN";
	if ($pthreads) {
		$pthreadsName = "Pthreads";
	} else {
		$pthreadsName = "NoPthreads";
	}
	return $pthreadsName;
}

sub getConcurrencyName()
{
	my $concurrencyName = "UNKNOWN";

	if ($mpi) {
		$concurrencyName="ConcurrencyMpi";
	} else {
		$concurrencyName="ConcurrencySerial";
	}
	return $concurrencyName;
}

sub getModelName()
{
	my $modelName = "UNKNOWN";

	if ($model=~/heisenberg/i) {
		$modelName="ModelHeisenberg";
	} elsif ($model=~/febasedsc/i) {
		$modelName = "ModelFeBasedSc";
	} elsif ($model=~/tjoneorbital/i) {
		$modelName = "TjOneOrbital";
	} elsif ($model=~/extendedhubbard1orb/i) {
		$modelName = "ExtendedHubbard1Orb";
	} elsif ($model=~/hubbard/i) {
		$modelName = "ModelHubbard"; # after extended
	}
	return $modelName;
}

sub getOperatorsName()
{
	my $operatorsName = "UNKNOWN";

	if ($model=~/extendedhubbard1orb/i) {
		$operatorsName = "OpsExtendedHubbard1Orb";
	} elsif ($model=~/hubbard/i) {
		$operatorsName = "OperatorsHubbard";
	} elsif ($model=~/heisenberg/i) {
		$operatorsName = "OperatorsHeisenberg";
	} elsif ($model=~/febasedsc/i) {
		$operatorsName="OperatorsFeAs";
	} elsif ($model=~/tjoneorbital/i) {
		$operatorsName="OperatorsTjOneOrbital";
	}
	return $operatorsName;
}

sub compilerName
{
	return "g++";
	my @tryThis = ("g++","g++4");
	my $ret;
	my $compiler;
	foreach my $comp (@tryThis) {
		my $ret = system("$comp &>2 /dev/null");
		if ($ret==0) {
			$compiler = $comp;
			last;
		} 
			
	}
	return $compiler if defined $compiler; 
	die "$0: No suitable compiler found\n";
}

