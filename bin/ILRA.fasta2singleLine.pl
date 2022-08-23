#!/usr/bin/env perl

# (c) Copyright: Genome Research Ltd. 2017
#
# Author: Thomas Dan Otto
# Licencse: GNU General Public License v3.0
#
# Description: put a multi fasta sequence into a single line

use strict;
use warnings;

my $largest = 0;
my $contig = '';


if (@ARGV != 1) {
	print "fasta2single.pl fasta > fasta.SL.fa\n\n" ;
	exit ;
}

my $filenameA = $ARGV[0];

open (IN, "$filenameA") or die "oops!\n" ;

	my $read_name = '' ;
	my $read_seq = '' ;

	while (<IN>) {
	    chomp ;

	    if (/^>(.+)/) {
		$read_name = $1 ;
		$read_seq = "" ;
		
		while (<IN>) {
		    chomp; 

			if (/^>(.+)/) {
			    
			    print ">$read_name\n$read_seq\n" ;


			    $read_name = $1 ;
			    $read_seq = "" ;



			}
			else {
			    chomp ;
			    s/\r//g;
			    $read_seq .= $_ ;

			}


		}

	    }
	}

close(IN) ;

print ">$read_name\n$read_seq\n" ;
