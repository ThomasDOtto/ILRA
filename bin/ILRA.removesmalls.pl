#!/usr/bin/env perl

## removesmalls.pl
## Unknown author
## Adapted by JLRuiz to be used in ILRA 2022, from https://www.biostars.org/p/79202/


use strict;
use warnings;

my $result = "";
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
        print ">$_" if($seqlen >= $minlen);
		$result .= $header."\n" if($seqlen <= $minlen);
    }
    local $/="\n";
}

open(my $fh, '>>', '../Excluded.contigs.fofn');
print $fh $result;
close $fh;

