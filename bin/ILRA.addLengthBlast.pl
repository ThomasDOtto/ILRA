#!/usr/bin/env perl
#
# File: IPA.addLengthBlast.pl
# Time-stamp: <13-Mar-2014 14:26:49 tdo>
# $Id: $
#
# Copyright (C) 2008 by Genome Research LTD
#
# Author: Thomas Dan Otto
# Licencse: GNU General Public License v3.0
# Description: put the length of fasta / fastq to a comparison output
# 
use strict;
use warnings;
use Data::Dumper;


my $fastq=shift;
my $reference=shift;

my $cigar = shift;

if ((!defined($cigar))) {
  die "usage: query cfasta <subject fasta> <reslut>\n";
  
}
my $ref_query=getFastaLength($fastq);
my $ref_subject=getFastaLength($reference);

#print Dumper $fastq;


open (F, $cigar) or die "Please give me cigar/ssaha...\n";

my @ar=<F>;


close(F);

my $res;


foreach (@ar) {
  chomp;
  
  my @a=split(/\s+/);
#  print "$a[1] $a[5] \n";
#  print " $$ref_query{$a[1]}\t$$ref_subject{$a[5]}    \n";
  
  $res .= $_."\t$$ref_query{$a[0]}\t$$ref_subject{$a[1]}\n";
  
}

open(F, "> $cigar.length") or die "Prob give me cigar/ssaha...\n";
print F $res;
close(F);


sub getFastqLength{
  my $name = shift;

  open (F, $name) or die "probl...\n";
  my @ar=<F>;
  
  close(F);

  my %h;
 
  
  for (my $i=0; $i<scalar(@ar);$i +=4 ) {
	chomp($ar[$i]);
	($name) = $ar[$i] =~ /^@(.*)$/;
	
	$h{$name}=length($ar[($i+1)]);
  }
  return \%h;
}
sub getFastaLength{
  my $filename = shift;

  open (F, $filename) or die "probl...\n";
  my @ar=<F>;
  
  close(F);
  my $Length=0;
  
  my %h;
  my $Counter=0;

  my $name;
  
  
  if ((0)) {
  
	for (my $i=0; $i<scalar(@ar);$i +=2 ) {
	  chomp($ar[$i]);
	  
	  ($name) = $ar[$i] =~ /^>(\S+)/;
	chomp($ar[($i+1)]);
#	print $1;
	
	  $h{$name}=length($ar[($i+1)]);
	}
  }
  
  else {
	foreach (@ar) {
	  if (/^>(\S+)/) {
		
		if ($Counter){ # not first
		  $h{$name}=$Length;
		}
		$name=$1;
#		print $name."\n";
		
		$Length =0;
		$Counter++;
	  }
	  else {
		chomp;
		$Length += length($_);
	  }
	  
	  
	  if ($Counter){
		$Counter++;
		$h{$name}=$Length;
	  }
	  
	}
  }
  return \%h;
}
