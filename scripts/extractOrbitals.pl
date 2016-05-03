#!/usr/bin/perl

use strict;
use warnings;

my ($orb1, $orb2, $orbitals) = @ARGV;

defined($orbitals) or die "USAGE: $0 orb1 orb2 orbitals\n";

my $total = ($orb1 == $orb2) ? 2 : 4;

print "#CONTINUEDFRACTIONCOLLECTION=$total\n";

my $offset = findOffset($orb1, $orb2, $orbitals);

print STDERR "$0: Called with $orb1 $orb2 $orbitals\n";
print STDERR "$0: Offset= $offset    total= $total\n";



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

sub findOffset
{
	my ($orb1, $orb2, $orbitals) = @_;
	($orb1 < $orbitals) or die "$0: $orb1 >= $orbitals\n";
	($orb2 < $orbitals) or die "$0: $orb2 >= $orbitals\n";
	my $offset = 0;
	
	for (my $x = 0; $x < $orbitals; ++$x) {
		for (my $y = $x; $y < $orbitals; ++$y) {
			return $offset if ($x == $orb1 and $y == $orb2);
			my $thisOffset = ($x == $y) ? 2 : 4;
			$offset += $thisOffset;
		}
	}

	return $offset;
}
