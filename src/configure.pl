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

#askQuestions();
my $model="hubbard";
my $modelLocation="Models/HubbardOneOrbital";

createMakefile();

createDriver();

#createObserverDriver();

sub createMakefile
{
	system("cp Makefile Makefile.bak") if (-r "Makefile");
	my $compiler = compilerName();
	my $headerFiles = join(' ', glob("Engine/*.h Models/*/*.h Geometries/*.h"));
	my @litProgFiles = glob("Engine/*.w Models/*/*.w Geometries/*.w");
	if ($modelLocation=~/extendedhubbard1orb/i) {
		$modelLocation = $modelLocation." -IModels/HubbardOneBand ";
	}
	my $litProgTargets = getLitProgTargets(\@litProgFiles);
	open(FOUT,">Makefile") or die "Cannot open Makefile for writing: $!\n";
print FOUT<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl
# Lanczos++ ($brand) by G.A.
# Platform: $platform
# MPI: $mpi

LDFLAGS =    $lapack  $gslLibs $pthreadsLib
CPPFLAGS = -Werror -Wall -I../PartialPsimag -IEngine -I$modelLocation -IGeometries -I$PsimagLite
EOF
if ($mpi) {
	print FOUT "CXX = mpicxx -O2 -DNDEBUG \n";
} else {
	print FOUT "CXX = $compiler -pg -O2 -DNDEBUG\n";
}
print FOUT<<EOF;
all: \$(EXENAME)
HEADERSH = $headerFiles

all: lanczos

lanczos: \$(HEADERSH)
	\$(CXX) -o lanczos \$(CPPFLAGS) lanczos.cpp \$(LDFLAGS)

clean:
	rm -f core* \$(EXENAME) *.o *.ii *.tt

$litProgTargets
######## End of Makefile ########

EOF
	close(FOUT);
	print STDERR "File Makefile has been written\n";
}

sub getLitProgTargets
{
	my ($array)=@_;
	my $x = "";
	my $litProgTool = "nuweb.pl -v -l  -s  -d ";
	foreach my $f (@$array) {
		my $fh = $f;
		$fh =~ s/\.w$/\.h/;
		$x = $x."$fh: $f\n";
		my $dir = $f;
		$dir =~ s/\/[^\/]+$/\//;
		my $fnd = $f;
		$fnd =~ s/$dir//;
		my $dirChange = computeBackwardMovements($dir);
		$x = $x."\t cd $dir; $dirChange$litProgTool $fnd\n";
		$x = $x."\n";
	}
	return $x;
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
	my $parametersName = "ParametersModelHubbard"; #getParametersName();
	my $pthreadsName = "NoPthreads.h"; #getPthreadsName();
	my $modelName = "HubbardLanczos"; #getModelName();
	my $operatorsName = getOperatorsName();
	
print FOUT<<EOF;
/* DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
 * This driver program was written by configure.pl
 * Lanczos++ ($brand) by G.A.*/
#include "$concurrencyName.h"
#include "ContinuedFraction.h"
#include "$modelName.h"
#include "$parametersName.h"
#include "Geometry.h"
#include "IoSimple.h"

using namespace Dmrg;
typedef double RealType;
typedef std::complex<RealType> ComplexType;
typedef $concurrencyName<RealType> ConcurrencyType;
typedef Geometry<RealType> GeometryType;
typedef $parametersName<RealType> ParametersModelType;
typedef IoSimple::In IoInputType;
typedef $modelName<RealType,ParametersModelType,GeometryType> ModelType;
typedef ContinuedFraction<ModelType,ConcurrencyType> ContinuedFractionType;
typedef ContinuedFractionType::TridiagonalMatrixType TridiagonalMatrixType;


ComplexType greenFunction(const TridiagonalMatrixType& ab,ComplexType z,RealType isign,
		const ContinuedFractionType& cf)
{
	//return ab.computeContinuedFraction(z,isign);
	return cf.continuedFraction(z,ab);
}

ComplexType greenFunction(RealType omega,RealType Eg,
		const TridiagonalMatrixType& abPlus,RealType normaPlus,
		const TridiagonalMatrixType& abMinus,RealType normaMinus,
		const ContinuedFractionType& cf)
{
	RealType isign = 1.0;
	RealType eps = 0.05;
	ComplexType z(omega + Eg,eps);
	ComplexType x = normaPlus * greenFunction(abPlus,z,isign,cf);
	RealType normalizer = 2.0;
	if (normaMinus!=0.0) {
		x -= normaMinus * greenFunction(abMinus,z,isign,cf);
		normalizer = 4.0;
	}
	return x/normalizer;
}

int main(int argc,char *argv[])
{
	//! setup distributed parallelization
	ConcurrencyType concurrency(argc,argv);

	//Setup the Geometry
	IoInputType io(argv[1]);
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
	size_t i = atoi(argv[2]);
	size_t j = atoi(argv[3]);

	ContinuedFractionType cf(model);

	//! get the g.s.:
	RealType Eg = cf.gsEnergy();
	std::cout<<"Energy="<<Eg<<"\\n";

	std::cout<<"gf(i="<<i<<",j="<<j<<")\\n";
	RealType normaPlus=0,normaMinus=0;
	TridiagonalMatrixType abPlus,abMinus;
	cf.getGreenFunction(abPlus,normaPlus,i,j,ContinuedFractionType::PLUS);
	if (i!=j) cf.getGreenFunction(abMinus,normaMinus,i,j,ContinuedFractionType::MINUS);
	for (int i=0;i<atoi(argv[4]);i++) {
		RealType omega = atof(argv[5])*i;
		std::cout<<omega<<" "<<std::imag(greenFunction(omega,Eg,abPlus,normaPlus,
				abMinus,normaMinus,cf))<<"\\n";
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
	my $parametersName="UNKNOWN";
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

