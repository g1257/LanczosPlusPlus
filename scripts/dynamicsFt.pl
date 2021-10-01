#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Math::Trig;
use Getopt::Long qw(:config no_ignore_case);

my $usage = "USAGE: $0 -r rootInputName  -n totalFiles -b begin -e end -s step -d delta [-A]\n";
my ($root, $total, $begin, $end, $delta, $step);
my $allPairs = 0; 
GetOptions('r=s' => \$root,
           'n=i' => \$total,
           'b=f' => \$begin,
           'e=f' => \$end,
           'd=f' => \$delta,
           's=f' => \$step,
           'A' => \$allPairs) or die "$usage\n";
my $myboolean = defined($root) && defined($total) && defined($begin) &&
defined($end) && defined($step) && defined($delta);

($myboolean) or die "$usage\n";

my @data;
for (my $m = 0; $m < $total; ++$m) {
	my $outName = doOneKmomentum($m, $root, $total, $allPairs);
	my $plotName = $root.$m.".dat";
	system("../../PsimagLite/drivers/continuedFractionCollection -d $delta -s $step -b $begin -e $end -f $outName > $plotName");
	print STDERR "$0: Written $plotName\n";
	my @a = loadPlot($plotName);
	$data[$m] = \@a;
}

printAll(\@data);

sub printAll
{
	my ($data) = @_;
	my $ptr0 = $data->[0];
	my @array0 = @$ptr0;
	my $omegas = scalar(@array0);
	for (my $i = 0; $i < $omegas; ++$i) {
		my $ptr2 = $array0[$i];
		my $omega2 = $ptr2->{"omega"};
		for (my $m = 0; $m < $total; ++$m) {
			my $k = findMomentumK($m, $total);
			my $ptr = $data->[$m];
			my @array = @$ptr;
			my $ptr3 = $array[$i];
			my $omega =  $ptr3->{"omega"};
			($omega == $omega2) or die "$0: Omegas don't match\n";
			my $value = $ptr3->{"data"};
			print "$k $omega $value->[0]\n";
		}

		print "\n";
	}
}

sub loadPlot
{
	my ($f) = @_;
	my @a;
	open(FILE, "<", $f) or die "$0: Cannot open $f : $!\n";
	while (<FILE>) {
		next if (/^#/);
		my @temp = split;
		(scalar(@temp) == 3) or die "$0: Expecting 3 values\n";
		my ($omega, $val0, $val1) = @temp;
		my $h = {"omega" => $omega, "data" => [$val0, $val1]};
		push @a, $h;
	}

	close(FILE);
	return @a;
}

sub doOneKmomentum
{
	my ($m, $root, $total, $allPairs) = @_;
	die "$0: All pairs option not implemented yet\n" if ($allPairs);
	my $sum = 0;
	my $text = "";
	my $total2 = ($allPairs) ? $total*$total : $total;
	for (my $i = 0; $i < $total2; ++$i) {
		my $finName = "$root$i.comb";
		my $n = procFile(\$text, $finName, $m, $i, $total);
		$sum += $n;
	}

	$text = "#CONTINUEDFRACTIONCOLLECTION=$sum\n$text\n";

	my $outName = $root."MomentumK".$m.".comb";
	open(FOUT, ">", "$outName") or die "$0: Cannot write to $outName : $!\n";
	print FOUT $text;
	close(FOUT);
	print STDERR "$0: Written to $outName $sum continued-fractions\n";
	return $outName;
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
			my $kmoment = findMomentumK($m, $total);
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

#Only for chains
sub findMomentumK
{
	my ($m, $total) = @_;
	return 2*pi*$m/$total;
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
