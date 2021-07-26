#!/usr/bin/perl -w

## removesmalls.pl
## Unknown author
## Adapted by JLRuiz from https://www.biostars.org/p/79202/
## 2021

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

open(my $fh, '>>', 'Excluded.contigs.fofn');
print $fh $result;
close $fh;

