#!/usr/bin/env perl
# Copyright (c) 2011-2015 Genome Research Ltd.
# Author: Thomas D. Otto <tdo@sanger.ac.uk>
#
# This file is part of ABACAS2.
#
# ABACAS2 is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;

my %h;

my $count=0;
while (<STDIN>){
	if (/^>/){
		if ($count){
			print "\n";
		}
		print ;
		$count++;
		}
	else {
		chomp;
		print ;
	}

}
