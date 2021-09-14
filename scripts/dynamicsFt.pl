#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Math::Trig;

my ($root, $total, $outRoot) = @ARGV;
defined($outRoot) or die "$0: rootInputName  totalFiles rootOutputName\n";

for (my $m = 0; $m < $total; ++$m) {
	doOneKmomentum($m, $root, $total, $outRoot);
}

sub doOneKmomentum
{
	my ($m, $root, $total, $outRoot) = @_;
	my $sum = 0;
	my $text = "";
	for (my $i = 0; $i < $total; ++$i) {
		my $finName = "$root$i.comb";
		my $n = procFile(\$text, $finName, $m, $i, $total);
		$sum += $n;
	}

	$text = "#CONTINUEDFRACTIONCOLLECTION=$sum\n$text\n";

	my $outName = "outRoot$m.comb";
	open(FOUT, ">", "$outName") or die "$0: Cannot write to $outName : $!\n";
	print FOUT $text;
	close(FOUT);
	print STDERR "$0: Written to $outName $sum continued-fractions\n";
}

sub procFile
{
	my ($text, $file, $m, $site, $total) = @_;
	my @factor;
	my $n;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		chomp;
		if (/CONTINUEDFRACTIONCOLLECTION=(.+$)/) {
			$n = $1;
			next;
		}

		if (/TSPCenter=(.+$)/) {
			my $center = $1;
			my $kmoment = 2*pi*$m/$total;
			my $arg = ($site - $center)*$kmoment;
			@factor = (cos($arg), sin($arg));
		}

		if (/#CFWeight=(.+$)/) {
			(scalar(@factor) == 2) or die "$0: TSPCenter not found in $file\n";
			my $weight = complexRealToReal($1);
			my @vals = ($factor[0]*$weight, $factor[1]*$weight);
			$$text .= "#CFWeight=(".$vals[0].",".$vals[1].")\n";
			next;
		}

		$$text .= "$_\n";
	}

	close(FILE);
	defined($n) or die "$0: CONTINUEDFRACTIONCOLLECTION not found in $file\n";
	return $n;
}

sub complexRealToReal
{
	my ($t) = @_;
	return $t if ($t =~ /^[\d\.eE\+\-]+$/);
	($t =~ /, *0\)$/) or die "$0: Expected imag part of weight to be zero not $t\n";
	$t =~ s/^\(//;
	$t =~ s/, *0\)$//;
	return $t;
}
