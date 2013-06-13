#!/usr/bin/perl 
=pod
// BEGIN LICENSE BLOCK
Copyright (c) 2009-2012, UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1.0.0]

---------------------------------------------------------------
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. 

Please see full open source license included in file LICENSE.
---------------------------------------------------------------
// END LICENSE BLOCK

=cut
use warnings;
use strict;

my $hasGsl = "no"; # say "no" here to remove GSL dependence

my $mpi=0;
my $platform="linux";
my $lapack="-llapack";
my $PsimagLite="../../PsimagLite";
my ($pthreads,$pthreadsLib)=(0,"");
my $brand= "v1.0";


my $gslLibs = " -lgsl  -lgslcblas ";
$gslLibs =" " if ($hasGsl=~/n/i);

#guessPlatform();

#welcome();

#my $model="";
#my $modelLocation="";
#my $stored="";
#askQuestions();

createMakefile();

#createDriver();


sub createMakefile
{
	system("cp Makefile Makefile.bak") if (-r "Makefile");
	my $compiler = compilerName();
	open(FOUT,">Makefile") or die "Cannot open Makefile for writing: $!\n";
	my $usePthreadsOrNot = " ";
	$usePthreadsOrNot = " -DUSE_PTHREADS " if ($pthreads);

print FOUT<<EOF;
# DO NOT EDIT!!! Changes will be lost. Modify configure.pl instead
# This Makefile was written by configure.pl
# Lanczos++ ($brand) by G.A.
# Platform: $platform
# MPI: $mpi

LDFLAGS =    $lapack  $gslLibs $pthreadsLib
CPPFLAGS = -Werror -Wall  -IEngine -IModels/Tj1Orb -IModels/Immm -IModels/HubbardOneOrbital -IModels/FeBasedSc -I$PsimagLite -I$PsimagLite/src $usePthreadsOrNot
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
	rm -f core* \$(EXENAME) *.o Makefile.dep

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

