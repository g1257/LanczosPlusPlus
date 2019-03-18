#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;

my ($template,$rootInput,$obs,$wbegin,$wend,$wstep,$wdelta,$orb1,$orb2,$orbitals,$spin) = @ARGV;
my $usage = "USAGE: $0 templateInput rootInput observable begin end step delta ";
$usage .= "orb1 orb2 orbitals [spin]";
defined($orbitals) or die "$usage\n";
defined($spin) or $spin = 0;
my $spins = "\"$spin,$spin\"";

my $total = readLabel($template,"TotalNumberOfSites");
my $centralSite = int($total/2) - 1;
my @data;

my $input = createInput($total);
system("./lanczos -f $input -g $obs -s $spins") if ($obs ne "0");

for (my $i = 0; $i < $total; ++$i) {
	system("perl ../scripts/extractOrbitals.pl $orb1 $orb2 $orbitals < $rootInput$i.comb > $rootInput$i.comb2");
	print STDERR "$0: Created $rootInput$i.comb2\n";
	system("echo \"#SITES $centralSite $i\" > $rootInput$i.cf");
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
my @intensities;
for (my $wi =0; $wi < $omegasTotal; ++$wi) {
	#print "$omegas[$wi]->[0] ";
	my @array;
	for (my $m = 0; $m < $total; ++$m) {
		my ($sumr, $sumi) = (0,0);
		my $q = 2.0*pi*$m/$total;
		for (my $i = 0; $i < $total; ++$i) {
			my $factor = ($i == $centralSite) ? 0.5 : 1;
			$_ = $data[$i];
			my @t = @$_;
			my $t = $t[$wi];
			my @temp = @$t;
			next unless ($temp[0] == $omegas[$wi]->[0]);
			my $arg = $q*($i-$centralSite);
			my $cosval = cos($arg);
			$sumi += $temp[1]*$cosval*$factor;
			$sumr += $temp[2]*$cosval*$factor;
		}

		#print "$sumi ";
		$array[$m] = {"kx" => $q, "value" => $sumi};
	}

	#print "\n";
	$intensities[$wi] = \@array;
}

printPgfPlot("$rootInput.pgfplots", \@intensities, \@omegas);


sub readData
{
	my ($a,$file) = @_;
	open(FILE, "<", "$file") or die "$0: Cannot open $file\n";
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
	my ($total) = @_;
	my $input = "$rootInput.inp";

	open(FILE, "<", $template) or die "$0: Cannot open $template : $!\n";
	open(FOUT, ">", "$input") or die "$0: Cannot write to $input : $!\n";
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

	my $c = int($total/2) - 1;
	for (my $i = 0; $i < $total; ++$i) {
		print FOUT "TSPSites 2   $c $i\n";
	}

	close(FOUT);

	print STDERR "$0: Created $input\n";
	return $input;
}

sub printPgfPlot
{
	my ($fout, $intensities, $omegas) = @_;
	open(FOUT, ">", "$fout") or die "$0: Cannot open $fout for writing: $!\n";
	# kx omega intensity
	# varying kx first
	for (my $wi =0; $wi < $omegasTotal; ++$wi) {
		my $omega = $omegas->[$wi]->[0];
		my $array = $intensities->[$wi];
		for (my $m = 0; $m < $total; ++$m) {
			my $h = $array->[$m];
			my $kx = $h->{"kx"};
			my $value = $h->{"value"};
			print FOUT "$kx $omega $value\n";
		}

		print FOUT "\n";
	}

	close(FOUT);
	print STDERR "$0: Written $fout\n";
}

sub readLabel
{
	my ($file, $label)=@_;
	open(FILE, "<", "$file") or die "$0: Cannot open $file\n";
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

