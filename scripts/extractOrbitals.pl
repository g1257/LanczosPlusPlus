#!/usr/bin/perl

use strict;
use warnings;

my ($orb1, $orb2, $orbitals) = @ARGV;
defined($orbitals) or die "USAGE: $0 orb1 orb2 orbitals\n";

my @keys;
while (<STDIN>) {
	if (/#INDEXTOCF /) {
		@keys = split;
		last;
	}
}

my ($offset,$total) = findOffsetAndTotal(\@keys,$orb1, $orb2, $orbitals);

print STDERR "$0: Called with $orb1 $orb2 $orbitals\n";
print STDERR "$0: Offset= $offset    total= $total\n";

print "#CONTINUEDFRACTIONCOLLECTION=$total\n";

while (<STDIN>) {
	last if (/#CONTINUEDFRACTIONCOLLECTION=/);
}

my $counter = 0;

while ($counter <= $offset) {
	$_ = <STDIN>;
	defined($_) or last;
	if (/^#Avector/) {
		$counter++;
	}
}

print "#Avector\n";

$counter = 0;
while ($counter <= $total) {
	$_ = <STDIN>;
	defined($_) or last;
	if (/^#Avector/) {
		$counter++;
	}

	print;
}

sub findOffsetAndTotal
{
	my ($keys,$orb1, $orb2, $orbitals) = @_;
	($orb1 < $orbitals) or die "$0: $orb1 >= $orbitals\n";
	($orb2 < $orbitals) or die "$0: $orb2 >= $orbitals\n";
	my ($offset,$total) = (0,0);
    	my $tempOff; 
	my $n = scalar(@$keys);
	if ($n == 0) {
		die "$0: Label INDEXTOCF not found in STDIN\n";
	}

	for (my $i = 1; $i < $n; ++$i) {
		$_ = $keys->[$i];
		my @temp = split/,/;
   
		my $b1 = ($orb1 == $temp[2] && $orb2 == $temp[3]);
		my $b2 = ($orb1 == $temp[3] && $orb2 == $temp[2]);
		if ($b1 || $b2) {
			$total++;
			$tempOff=$i;
		} 
	}

	if ($tempOff < $total) {
		die "$0: Internal error $tempOff < $total\n";
	}

	$offset = $tempOff - $total;
	return ($offset,$total);
}
