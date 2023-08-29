#!/usr/bin/env perl

## removesmalls.pl
## Unknown author
## Adapted by JLRuiz to be used in ILRA 2022, from https://www.biostars.org/p/79202/


use strict;
use warnings;

my $result = "";
my $result2 = "";
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
        print ">$_" if($seqlen >= $minlen);
		$result .= $header."\n" if($seqlen <= $minlen);
    }
    local $/="\n";
}

open(my $fh, '>>', '../Excluded.contigs.fofn');
print $fh $result;
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
print "\nR code in case you want to get the distribution:";
my $message = q|Rscript -e 'library(ggplot2); data <- read.table("Contigs_length.txt", header=FALSE, stringsAsFactors=FALSE); colnames(data) <- c("Contig", "Value")
Q1 <- quantile(data$Value, 0.25); Q3 <- quantile(data$Value, 0.75); data$Value_Kbp <- data$Value/1000
p <- ggplot(data, aes(x = "", y = Value_Kbp)) + geom_violin(fill = "lightblue") + geom_boxplot(width = 0.1, fill = "white", color = "black") +
annotate("text", x = 0.2, y = Q1/1000, label = paste("Q1: ", round(Q1, 2), " bp"), hjust=0) + annotate("text", x = 0.2, y = Q3/1000, label = paste("Q3: ", round(Q3, 2), " bp"), hjust=0) +
theme_classic() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
pdf("Contig_length_distribution.pdf");print(p);dev.off()'|;
print $message;
close $fh3;


