#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Math::Trig;

my ($root, $total) = @ARGV;
defined($total) or die "$0: rootInputName  totalFiles\n";

my $sum = 0;
my $text = "";
for (my $i = 0; $i < $total; ++$i) {
	my $arg = 2*pi*$i/$total;
	my @factor = (cos($arg), sin($arg));
	my $finName = "$root$i.comb";
	my $n = procFile(\$text, $finName, \@factor);
	$sum += $n;
}

$text = "#CONTINUEDFRACTIONCOLLECTION=$sum\n$text\n";

print $text;

sub procFile
{
	my ($text, $file, $factor) = @_;
	my $n;
	open(FILE, "<", $file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		chomp;
		if (/CONTINUEDFRACTIONCOLLECTION=(.+$)/) {
			$n = $1;
			next;
		}

		if (/#CFWeight=(.+$)/) {
			my $weight = $1;
			my @vals = ($factor->[0]*$weight, $factor->[1]*$weight);
			$$text .= "#CFWeight=($vals[0], $vals[1])\n";
			next;
		}
		
		$$text .= "$_\n";
	}
	
	close(FILE);
	defined($n) or die "$0: CONTINUEDFRACTIONCOLLECTION not found in $file\n";
	return $n;
}
