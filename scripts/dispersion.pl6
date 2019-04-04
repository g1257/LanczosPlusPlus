#!/usr/bin/perl6

use v6;

sub MAIN(Int $L, $t, $r, Str $bc = "open")
{
	my $tpr = $t + $r;
	my $tmr = $t - $r;
	my @disp;
	for 0..^$L -> Int $mm {
		my $m = ($bc eq "periodic") ?? $mm !! $mm + 1;
		my $sk = -2.0*cos(findKdependence($m, $L, $bc));
		my $e = $tpr*$sk;
		push @disp, $e;
		$e = $tmr*$sk;
		push @disp, $e;
	}

	my @disp2 = @disp.sort;
	say "Eigenvalues: ";
	say @disp2.join(' ');
	say "Full Energies: ";
	say sumEnergies(@disp2).join("\n");
}

sub sumEnergies(@disp)
{
	my @a;
	my $sum = 0;
	my $n = @disp.elems;
	for 0..^$n -> Int $ind {
		$sum += @disp[$ind];
		my Int $jnd = $ind + 1;
		push @a, "$jnd\t$sum";
	}

	return @a;
}

sub findKdependence($m, $L, $bc)
{
	return ($bc eq "periodic") ?? kperiodic($m, $L) !! kopen($m, $L);
}

sub kperiodic($m, $L)
{
	2.0*$m*pi/$L;
}

sub kopen($m, $L)
{
	$m*pi/($L + 1);
}


