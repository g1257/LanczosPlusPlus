#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;

my ($template,$rootInput,$obs,$wbegin,$wend,$wstep,$wdelta) = @ARGV;
my $usage = "USAGE: $0 templateInput rootInput observable begin end step delta";
defined($wdelta) or die "$usage\n";

my $total = readLabel($template,"TotalNumberOfSites");
my $centralSite = int($total/2) - 1;
my @data;

for (my $i = 0; $i < $total; ++$i) {
	my $input = createInput($i);
	system("./lanczos -f $input -g $obs &> $rootInput$i.comb");
	print STDERR "$0: Created $rootInput$i.comb\n";
	system("echo \"#SITES $centralSite $i\" > $rootInput$i.cf");
	my $cmd = "../../PsimagLite/drivers/continuedFractionCollection ";
	$cmd .= "-f $rootInput$i.comb -b $wbegin ";
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
	for (my $m = 0; $m < $total; ++$m) {
		my ($sumr, $sumi) = (0,0);
		my $q = 2.0*pi*$m/$total;
		for (my $i = 0; $i < $total; ++$i) {
			$_ = $data[$i];
			my @t = @$_;
			my $t = $t[$wi];
			my @temp = @$t;
			next unless ($temp[0] == $omegas[$wi]->[0]);
			my $arg = $q*($i-$centralSite);
			my $cosval = cos($arg);
			my $sinval = sin($arg);
			$sumr += $temp[1]*$cosval - $temp[2]*$sinval;
			$sumi += $temp[1]*$sinval + $temp[2]*$cosval;
		}

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
	my $site1 = $centralSite;
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
