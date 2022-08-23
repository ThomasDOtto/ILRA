#!/usr/bin/env perl
#
# File: IPA.deleteEntryinBlast.pl
# Time-stamp: <25-Feb-2013 17:19:57 tdo>
# $Id: $
#
# (c) Copyright: Genome Research Ltd. 2017
#
# Author: Thomas Dan Otto
# Licencse: GNU General Public License v3.0
## Copyright (C) 2013 by Pathogene Group, Sanger Center
#
# Description: Deletes a list of query from a blast file
#
use strict;
use warnings;

my $f=shift;


my $ref_l=getList($f);


while (<STDIN>) {
  my @ar=split(/\t/);

  if ( (!defined($$ref_l{$ar[0]})) && (! defined($$ref_l{$ar[1]}))) {
	print ;
  }
}

sub getList{
  open F, shift or die "gimme list\n";
  my %h;
 my %s;
  while (<F>) {
	chomp;
	$h{$_}=1;
#	my @ar=split(/\t/);
#	$h{$ar[0]}=$ar[1];
  }
  return \%h
}
