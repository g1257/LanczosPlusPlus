#!/usr/bin/perl
=pod
Copyright (c) 2009-2015, UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

=cut
use warnings;
use strict;

use Getopt::Long qw(:config no_ignore_case);
use lib "../../PsimagLite/scripts";
use NewMake;
use lib ".";
use PsiTag;

my ($flavor, $lto) = (NewMake::noFlavor() , 0);
my $usage = "USAGE: $0 [-f flavor] [-lto] [-c config]\n";
my $config;

GetOptions('f=s' => \$flavor,
           'lto' => \$lto,
           'c=s' => \$config) or die "$usage\n";

my $gccdash = "";
if ($lto == 1) {
	$gccdash = "gcc-";
	$lto = "-flto";
} else {
	$lto = "";
}

my @configFiles = ("../TestSuite/inputs/ConfigBase.psiTag");
push @configFiles, $config if (defined($config));

my %thermalD = (name => "thermal", dotos => "thermal.o");
my %lorentzianD = (name => "lorentzian", dotos => "lorentzian.o");
my %ld0 = (name => "LanczosDriver0", aux => 1);
my %ld1 = (name => "LanczosDriver1", aux => 1);
my %ld2 = (name => "LanczosDriver2", aux => 1);
my %ld3 = (name => "LanczosDriver3", aux => 1);
my %ld4 = (name => "LanczosDriver4", aux => 1);
my %ld5 = (name => "LanczosDriver5", aux => 1);
my $dotos = "lanczos.o LanczosDriver0.o LanczosDriver1.o LanczosDriver2.o";
$dotos .= " LanczosDriver3.o LanczosDriver4.o LanczosDriver5.o";
my %lanczosD = (name => "lanczos", dotos => $dotos);

my @drivers = (\%thermalD,
	\%lorentzianD,
	\%ld0,
	\%ld1,
	\%ld2,
	\%ld3,
	\%ld4,
	\%ld5,
	\%lanczosD);

my %args;
$args{"CPPFLAGS"} = $lto;
$args{"LDFLAGS"} = $lto;
$args{"flavor"} = $flavor;
$args{"code"} = "Lanczos++";
$args{"configFiles"} = \@configFiles;

createMakefile(\@drivers, \%args);

sub createMakefile
{
	my ($drivers, $args) = @_;
	NewMake::backupMakefile();

	my $fh;
	open($fh, ">", "Makefile") or die "Cannot open Makefile for writing: $!\n";

	NewMake::main($fh, $args, \@drivers);

	close($fh);
	print STDERR "$0: File Makefile has been written\n";
}

