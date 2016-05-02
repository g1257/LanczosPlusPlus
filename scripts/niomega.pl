#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;

my ($template,$rootInput,$wbegin,$wend,$wstep,$wdelta,$orb1,$orb2,$orbitals) = @ARGV;
my $obs = "c";
my $usage = "USAGE: $0 templateInput rootInput begin end step delta orb1 orb2 orbitals";
defined($orbitals) or die "$usage\n";

my $total = readLabel($template,"TotalNumberOfSites");
my @data;

for (my $i = 0; $i < $total; ++$i) {
	my $input = createInput($i);
	system("./lanczos -f $input -g $obs &> $rootInput$i.comb");
	print STDERR "$0: Created $rootInput$i.comb\n";
	system("perl ../scripts/extractOrbitals.pl $orb1 $orb2 $orbitals < $rootInput$i.comb > $rootInput$i.comb2");
	print STDERR "$0: Created $rootInput$i.comb2\n";
	system("echo \"#SITES $i $i\" > $rootInput$i.cf");
	my $cmd = "../../PsimagLite/drivers/continuedFractionCollection ";
	$cmd .= "-f $rootInput$i.comb2 -b $wbegin ";
	$cmd .= " -e $wend -s $wstep -d $wdelta >> $rootInput$i.cf";
	system($cmd);
	print STDERR "$0: Created $rootInput$i.cf\n";
	my @temp;
	readData(\@temp,"$rootInput$i.cf");
	$data[$i] = \@temp;
}

$_ = $data[0];
my @omegas = @$_;
my $omegasTotal = scalar(@omegas);
print STDERR "#Omegas=".$omegasTotal."\n";
for (my $wi =0; $wi < $omegasTotal; ++$wi) {
	print "$omegas[$wi]->[0] ";
	my ($sumr, $sumi) = (0,0);
	for (my $i = 0; $i < $total; ++$i) {
		$_ = $data[$i];
		my @t = @$_;
		my $t = $t[$wi];
		my @temp = @$t;
		next unless ($temp[0] == $omegas[$wi]->[0]);
		$sumi = $temp[1];

		print "$sumi ";
	}

	print "\n";
}


sub readData
{
	my ($a,$file) = @_;
	open(FILE,"$file") or die "$0: Cannot open $file\n";
	my $value;
	my $counter = 0;
	while(<FILE>) {
		chomp;
		next if (/^#/);
		my @temp = split;
		next if (scalar(@temp) != 3);
		$a->[$counter++] = \@temp;
	}

	close(FILE);

	print STDERR "$0: Read $counter rows from $file\n";
}

sub createInput
{
	my ($ind) = @_;
	my $input = "$rootInput$ind.inp";
	my $site1 = $ind;
	my $site2 = $ind;

	open(FILE,$template) or die "$0: Cannot open $template : $!\n";
	open(FOUT,">$input") or die "$0: Cannot write to $input : $!\n";
	while(<FILE>) {
		next if (/^#/);
		if (/\$([a-zA-Z0-9\[\]]+)/) {
				my $name = $1;
				my $str = "\$".$name;
				my $val = eval "$str";
				defined($val) or die "$0: Undefined substitution for $name\n";
				s/\$\Q$name/$val/g;
		}

		print FOUT;
	}

	close(FILE);
	close(FOUT);

	print STDERR "$0: Created $input\n";
	return $input;
}

sub readLabel
{
	my ($file, $label)=@_;
	open(FILE,"$file") or die "$0: Cannot open $file\n";
	my $value;
	while(<FILE>) {
		chomp;
		if (/^$label=(.*$)/) {
			$value=$1;
			last;
		}
	}

	close(FILE);

	defined($value) or die "$0: Label $label not found in $file\n";
	return $value;
}

