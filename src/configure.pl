#!/usr/bin/perl
=pod
Copyright (c) 2009-2014, UT-Battelle, LLC
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

=cut
use warnings;
use strict;

use lib "../../PsimagLite/scripts";
use Make;

my $mpi=0;
my $platform=guessPlatform();
my @drivers = ("lanczos");
my $lapack=Make::findLapack();
my $PsimagLite="../../PsimagLite";
my ($pthreads,$pthreadsLib)=(0,"");
my $brand= "v1.0";
my $build="production";
my $floating="";

system("make clean");

createMakefile();

sub guessPlatform
{
	my $platform="Linux";
	$platform="Darwin" if (isAMac());
	return $platform;
}

sub isAMac
{
	open(PIPE,"uname -a |grep -i Darwin") or return 0;
	$_=<PIPE>;
	close(PIPE);

	return 1 unless ($_ eq "" or $_ eq "\n");
	return 0;
}

sub createMakefile
{
	unlink("Engine/Version.h");
	system("cp Makefile Makefile.bak") if (-r "Makefile");
	my $compiler = "g++ ";
	$compiler = " mpicxx " if ($mpi);
	my $fh;
	open($fh,">Makefile") or die "Cannot open Makefile for writing: $!\n";
	my $usePthreadsOrNot = " ";
	$usePthreadsOrNot = " -DUSE_PTHREADS " if ($pthreads);
	my $optimizations = " -O3 -DNDEBUG ";
	$optimizations = " -g3 -D_GLIBCXX_DEBUG " if ($build eq "debug");
	$optimizations .= " -g3 " if ($build eq "callgrind");
	my $strip = "strip ";
	$strip = " true " if ($build eq "debug" or $build eq "callgrind");


	my $cppflags= "-Werror -Wall -Wstrict-overflow=5  -IEngine  ";
	$cppflags .= "  -I$PsimagLite/src -I$PsimagLite $usePthreadsOrNot $floating";

	Make::make($fh,\@drivers,"DMRG++",$platform,$mpi,"$lapack $pthreadsLib",
	"$compiler $optimizations",$cppflags,$strip,"Engine/Version.h",
	"Engine/Version.h gitrev","");
	local *FH = $fh;
print FH<<EOF;

Engine/Version.h: gitrev
	./gitrev > Engine/Version.h

EOF

	close($fh);
	print STDERR "$0: File Makefile has been written\n";
}


