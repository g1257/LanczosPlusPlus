#!/usr/bin/perl

use strict;
use warnings;

my ($output,$templateInput,$total) = @ARGV;
defined($total) or die "USAGE: $0 outputFile inputTemplate total\n";

if (-e "$output") {
	print "$output file exists, delete? ";
	$_=<STDIN>;
	chomp;
	$_ eq "y" or $_ eq "yes" or die "$0: Aborted by user\n";
}

unlink($output);
my $totalSector = "#TotalSectors=".($total+1);
system("echo \"$totalSector\" > $output");

for (my $spc = 0; $spc <= $total; ++$spc) {
	my $ndown = $total - $spc;
	my $input = createInput($spc,$total,$spc,$ndown);
	runThis($output,$input);
	print STDERR "$0: Appended run $spc in $output\n";
}

sub createInput
{
	my ($spc,$total,$nup,$ndown) = @_;
	my $inputFile = "Input$spc.inp";

	my $hubbardU = createVector($total,1);
	my $potentialV = createVector($total,0);
	open(FOUT,">$inputFile") or die "$0: Cannot write to $output: $!\n";
	open(FILE,"$templateInput") or die "$0: Cannot open $templateInput: $!\n";

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

