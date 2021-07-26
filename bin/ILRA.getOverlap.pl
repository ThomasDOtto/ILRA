#!/usr/bin/perl -w
#
# (c) Copyright: Genome Research Ltd. 2017
#
# Author: Thomas Dan Otto
# Licencse: GNU General Public License v3.0
#
#
#
# Author: Thomas Otto
#
# Description: Get the overlap from a m8 blast output.
#

use strict;


## it will get an m8 blast and return the overlap / hit to query


my %seen;
my %qLength;

while (<STDIN>) {
  chomp;
  
  my ($q,$s,$id,$ov,$one,$two,$qStart,$qEnd,$sStart,$sEnd,
	  $score,$value,$qLength,$sLength)=split(/\t/);
  $qLength{$q}=$qLength;
  
  for my $i ($qStart..$qEnd){
	$seen{$q}{$i}=1;      
	
  }
}

foreach my $g (keys %qLength) {
  my $overlap=0;
  
  foreach my $pos (keys %{ $seen{$g} }) {
	$overlap++
  }
  print "$g\t$qLength{$g}\t$overlap\t".int($overlap*100/$qLength{$g})."\n";
  
}


sub notSeen{
        my $q= shift;
        my $genome=shift;
        my $s=shift;
        my $e=shift;
        
        for ($s..$e){
                if (defined($seen{$q}{$genome}{$_})){
                        return 0        
                }       
        }       
        
        return 1;       
}
  
