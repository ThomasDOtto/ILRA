#!/usr/bin/env perl

## removesmalls.pl
## Unknown author
## Adapted by JLRuiz to be used in ILRA 2022, from https://www.biostars.org/p/79202/


use strict;
use warnings;

my $result = "";
my $result2 = "";
my $result3 = "";
my $minlen = shift or die "Error: `minlen` parameter not provided\n";
{
    local $/=">";
    while(<>) {
        chomp;
        next unless /\w/;
        s/>$//gs;
        my @chunk = split /\n/;
        my $header = shift @chunk;
        my $seqlen = length join "", @chunk;
        $result2 .= $header."\t".$seqlen."\n";
        $result3 .= $header."\t".$seqlen."\n" if($seqlen < $minlen);
        print ">$_" if($seqlen >= $minlen);
		$result .= $header."\n" if($seqlen <= $minlen);
    }
    local $/="\n";
}

open(my $fh, '>>', '../Excluded.contigs.fofn');
print $fh $result3;
close $fh;

open(my $fh2, '>>', 'Contigs_length.txt');
print $fh2 $result2;
close $fh2;

## Compute and print Q1 and Q3:
# Assuming $result2 is your data
my @lines = split("\n", $result2);
my @values;

foreach my $line (@lines) {
    my ($name, $value) = split("\t", $line);
    push @values, $value;
}

@values = sort { $a <=> $b } @values;

sub quantile_type7 {
    my ($data_ref, $prob) = @_;
    my @sorted_data = sort {$a <=> $b} @$data_ref;
    my $n = scalar @sorted_data;

    my $h = ($n - 1) * $prob + 1;
    my $j = int($h);
    my $g = $h - $j;

    if ($j == $n) {
        return $sorted_data[$j - 1];
    } else {
        return (1 - $g) * $sorted_data[$j - 1] + $g * $sorted_data[$j];
    }
}


my $q1 = quantile_type7(\@values, 0.25);
my $q3 = quantile_type7(\@values, 0.75);

## Print messages:
open(my $fh3, '>>', 'Contigs_length_stats.txt');
select $fh3;
print "Contigs length Q1: $q1\n";
print "Contigs length Q3: $q3\n";
close $fh3;


