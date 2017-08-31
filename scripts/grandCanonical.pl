#!/usr/bin/perl

use strict;
use warnings;

my ($model,$output,$templateInput,$ntotal) = @ARGV;
defined($templateInput) or die "USAGE: $0 model outputFile inputTemplate";

my $total = readLabel($templateInput,"TotalNumberOfSites=");
print STDERR "TotalNumberOfSites=$total\n";

if ($model eq "canonical") {
	defined($ntotal) or die "USAGE: $0 model outputFile inputTemplate total [ntotal]\n";
} else {
	$ntotal = $total;
}

if (-e "$output") {
	unlink($output);
}

my $sectors = 0;
for (my $nup = 0; $nup <= $total; ++$nup) {
	for (my $ndown = 0; $ndown <= $total; ++$ndown) {
		if ($model eq "Heisenberg" || $model eq "canonical") {
			next if ($nup + $ndown != $ntotal);
		}

		if ($model eq "tj") {
			next if ($nup + $ndown >= $total);
		}

		$sectors++;
	}
}

my $totalSector = "#TotalSectors=$sectors";
system("echo \"$totalSector\" > $output");

my $counter = 0;
for (my $nup = 0; $nup <= $total; ++$nup) {
	for (my $ndown = 0; $ndown <= $total; ++$ndown) {
		if ($model eq "Heisenberg"|| $model eq "canonical") {
			next if ($nup + $ndown != $ntotal);
		}

		if ($model eq "tj") {
			next if ($nup + $ndown >= $total);
		}

		$counter++;
		my $input = createInput($counter,$total,$nup,$ndown);
		runThis($output,$input);
		print STDERR "$0: Appended run $nup $ndown in $output\n";
	}
}

sub createInput
{
	my ($spc,$total,$nup,$ndown) = @_;
	my $inputFile = "Input$spc.inp";

	my $hubbardU = createVector($total,1);
	my $potentialV = createVector($total,0);
	open(FOUT, ">", "$inputFile") or die "$0: Cannot write to $output: $!\n";
	open(FILE, "<", "$templateInput") or die "$0: Cannot open $templateInput: $!\n";

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
	print STDERR "$0: File $inputFile has been written\n";
	return $inputFile;
}

sub runThis
{
	my ($output,$input) = @_;
	my $cmd = "./lanczos -f $input >> $output 2>/dev/null";
	print STDERR "$0: Execing $cmd\n";
	system($cmd);
}

sub createVector
{
	my ($n,$val) = @_;
	my $a = "$n ";
	for (my $i = 0; $i < $n; ++$i) {
		$a .= "$val ";
	}

	return $a;
}

sub readLabel
{
	my ($file,$label)=@_;
	my $ret;
	open(FILE, "<", $file) or die "Cannot open $file: $!\n";
	while(<FILE>) {
		chomp;
		if (/^$label(.*$)/) {
			$ret=$1;
			last;
		}
	}

	close(FILE);

	defined($ret) or die "readLabel: Not found $label in $file\n";
	return $ret;
}
