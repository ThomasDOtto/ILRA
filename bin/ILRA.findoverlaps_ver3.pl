#!/usr/bin/perl -w

##################################################
### INPUT:
### 1. file with blast results, length of sequences in last two columns
### 2. bam file with reads mapped to the assembly (preferably only proper pairs)
### 3. fasta file including the assembly
### 4. identifier for output files
###
### OUTPUT:
### 1. mergedseq_output.fasta: fasta file with the new assembly
### 2. result_output.txt: information which sequences have been merged
### 3. log_output.txt: number of results for the different filtering steps
### 4. contained_output.fasta: fasta file with contained sequences
### 5. notcovered_output.fasta: fasta file with sequences that have no coverage
#################################################

use strict; 
use Data::Dumper;

### get arguments

my $blastfile = shift;
my $bamfile  = shift;
my $fastafile = shift;
my $out = shift;


if (! defined($out)){
    print "usage: findoverlaps_ver3.pl <blastoutput> <bamfile> <fastafile> <output_name> \n";
    exit;
}


### names for output files
my $newfastafile = "mergedseq_".$out.".fasta";
my $resultfile = "results_".$out.".txt";
my $logfile = "log_".$out.".txt";
my $containedfasta = "contained_".$out.".fasta";
my $notcovered = "notcovered_".$out.".fasta";



open(LOG, ">$logfile") or die "problem writing logfile $logfile: $!";


my $ref_coverage = loadCoverage($bamfile);
my $ref_notcovered = notCovered($ref_coverage);

my $ref_contained = getcontained($blastfile, $ref_coverage, $ref_notcovered);

my $size = @$ref_notcovered;
print "$size\n";

my @ignore = (@$ref_notcovered, @$ref_contained);    #combine the two arrays with sequences that should be ignored

my $ref_blast = loadBlast($blastfile, \@ignore);
my $number = keys %$ref_blast;
print LOG "Blast results after filtering out contained and not covered sequences\nresults: $number\n\n";


my $ref_overlaps = findOverlaps($ref_blast);
my $ref_overlaps = evaluateOverlaps($ref_overlaps, $ref_coverage);
my $number = keys %$ref_overlaps;
print LOG "Second filtering step: evaluate the coverage\nresults: $number\n\n";

my $ref_overlaps = getmerging($ref_overlaps, $fastafile);


my $ref_seq = readseq($fastafile, $ref_overlaps, $ref_contained, $ref_notcovered, $newfastafile, $containedfasta, $notcovered);
my $ref_mergedseq = mergeoverlaps($ref_overlaps, $ref_seq, $resultfile);

writenewfasta($ref_mergedseq, $newfastafile);



###################################
#### subroutines
###################################



exit 0;


sub notCovered{
    my $coverage = shift;
    
    ###array to save names of sequences that have no coverage
    my @notcovered;
 
    ###get overall average coverage
    ###calculate the mean coverage of the genome
    my $sum = 0;
    my $counter = 1;
    LOOP: foreach my $key1 (keys %$coverage){          #loop over sequences
      foreach my $key2 (keys %{$$coverage{$key1}}){  #loop over position
	  if ($counter > 1000000){                   #only use first million bases for calculation
	      last LOOP;
	  }
	  $sum = $sum + $$coverage{$key1}{$key2};
	  $counter++;
      }

    }
    my $mean = $sum / $counter;
#    print "$mean";

    ###loop over sequences
    foreach my $key1 (keys %$coverage){
	my $sum = 0;
	my $counter = 1;
	###loop over position
	foreach my $key2 (keys %{$$coverage{$key1}}){
	    $sum = $sum + $$coverage{$key1}{$key2};
	    $counter++;
	}
	
	my $average = $sum / $counter;
	###check coverage
	if ($average <= $mean/10){
	    push(@notcovered, $key1);
#	    print "$key1\t$average\n";
	}
    }

    return \@notcovered;
}


sub getcontained{
    my $blastfile = shift;
    my $coverage = shift;
    my $notcovered = shift;
 
    my $align = loadBlast($blastfile, $notcovered);
   
    ###filter for contained sequences
    for my $key (keys %$align){
	###remove alignments where query and subject are the same
	if ($$align{$key}{queryID} eq $$align{$key}{subjectID}){
	    delete $$align{$key};
	    next;
	}

	my $dev = 50;          #allowed deviation

	###check if subject contained
	if ($$align{$key}{s_start} <= $dev && $$align{$key}{s_end} >= $$align{$key}{s_length} - $dev) {   #not reversed
	    $$align{$key}{contained} = "subject";
	} 
	elsif ($$align{$key}{s_end} <= $dev && $$align{$key}{s_start} >= $$align{$key}{s_length} - $dev){  #reversed
	    $$align{$key}{contained} = "subject";
	    $$align{$key}{reversed} = 2;
	}

	###check if query contained
	elsif ($$align{$key}{q_start} <= $dev && $$align{$key}{q_end} >= $$align{$key}{q_length} - $dev){  
	    $$align{$key}{contained} = "query";
	}
	elsif ($$align{$key}{q_end} <= $dev && $$align{$key}{q_start} >= $$align{$key}{q_length} - $dev){
	    $$align{$key}{contained} = "query";
	    $$align{$key}{reversed} = 1;
	}

	else {
	    delete $$align{$key};
	}
    }

    ###filter double results (results for same alignment but switched query/subject)
    foreach my $key1 (keys %$align){
	foreach my $key2 (keys %$align){
	    unless (exists $$align{$key1} && exists $$align{$key2}){  #check if keys still exists
		next;
	    }

	    my $delete = 0;

	    if($$align{$key1}{subjectID} eq $$align{$key2}{queryID} 
	       && $$align{$key1}{queryID} eq $$align{$key2}{subjectID}){
		if ($$align{$key1}{reversed} == 0 && $$align{$key2}{reversed} == 0){
		    if ($$align{$key1}{s_start} == $$align{$key2}{q_start}
			&& $$align{$key1}{s_end} == $$align{$key2}{q_end}
			&& $$align{$key1}{q_start} == $$align{$key2}{s_start}
			&& $$align{$key1}{q_end} == $$align{$key2}{s_end}){
			$delete = 1;
		    }
		}else{
		    if ($$align{$key1}{s_start} == $$align{$key2}{q_end}
			& $$align{$key1}{s_end} == $$align{$key2}{q_start}
			& $$align{$key1}{q_start} == $$align{$key2}{s_end}
		        & $$align{$key1}{q_end} == $$align{$key2}{s_start}){
			$delete = 1;
		    }
		}
	    }

	    if($delete == 1){
		delete $$align{$key2};
	    }
	}
    }

    my $number = keys %$align;
    print LOG "Contained before checking coverage: $number\n";

    ###check if coverage is lower than at the rest of the genome
    evaluateOverlaps($align, $coverage);

    my $number = keys %$align;
    print LOG "Contained after checking coverage: $number\n";

    my @contained;
    foreach my $key (keys %$align){
	if ($$align{$key}{contained} eq "subject"){
	    push(@contained, $$align{$key}{subjectID});
	}
	elsif ($$align{$key}{contained} eq "query"){
	    push(@contained, $$align{$key}{queryID});
	}
    }

    return \@contained;
}


sub loadBlast{
    
    my $filename = shift;
    my $ignore = shift;

    ###save ignored sequences as keys in a hash
    my %ignore = map { $_ => 1 } @$ignore;

    open (blast, $filename) or die "problem opening blastfile $filename : $!\n";

    ### hash to store the alignemnts
    my %align;

    ### read data and store alignments in hash
    my $i = 1;
    LINE: while (<blast>){
		unless(/^#/){		#ignore comments
                	my @values = split("\t");
			
			for (@values){          #remove whitespaces at the beginning/end
               			s/^\s+|\s+$//;
        		}

			#check if subject or query should be ignored, if yes: skip
			if (exists($ignore{$values[0]}) || exists($ignore{$values[1]})){  
			next LINE;    
			}

                	$align{$i}= {
                        	queryID => $values[0],
                       		subjectID => $values[1],
				#identity => $values[2],
				alignment_length => $values[3],
				#mismatches => $values[4],
				#gap_openings => $values[5],
				q_start => $values[6],
				q_end => $values[7],
				s_start => $values[8],
				s_end => $values[9],
				#e_value => $values[10],
				#bit_score => $values[11],
				q_length => $values[12],
				s_length => $values[13]
                	};
                	$i++;
        	}
    }
    close blast;

    return \%align;
}


sub findOverlaps{
	my $align = shift;

	###declaring what deviation from the start/end of the sequences should be allowed	
	my $dev = 50;


	###remove alignments where query and subject are the same
	foreach my $key (keys %$align){
	    if ($$align{$key}{queryID} eq $$align{$key}{subjectID}){
		delete $align->{$key};
		next;
	    }
	}


	###filter for overlap alignments
	foreach my $key (keys %$align){

	    my $delete = 1;

	    ### filter for the different cases of overlap alignments
	    ### if an alignment is found: save information about order and if one seq is reversed for later

	    #first case: alignment at the beginning of one sequence and at the end of the other
	    if ($$align{$key}{q_end} >= ($$align{$key}{q_length} - $dev) && $$align{$key}{s_start} <= (1 + $dev)){
		#alignment on the right of query and left on subject
		$delete = 0;
		$$align{$key}{q_align} = "right";
		$$align{$key}{s_align} = "left";
	    } elsif($$align{$key}{q_start} <= (1 + $dev) && $$align{$key}{s_end} >= ($$align{$key}{s_length} - $dev)){
		#alignment on the left of query and right on subject
		$delete = 0;
		$$align{$key}{q_align} = "left";
		$$align{$key}{s_align} = "right";
	    }

	    #second case: alignment at the beginning or end of both sequences
	    #(one alignment reversed)
	    
	    #query sequence reversed
	    if(($$align{$key}{q_start} >= $$align{$key}{q_end})){
		if ($$align{$key}{q_start} >= ($$align{$key}{q_length} - $dev) && $$align{$key}{s_end} >= ($$align{$key}{s_length} - $dev)){
		    $delete = 0;
		    $$align{$key}{reversed} = 1;
		    $$align{$key}{q_align} = "right";
		    $$align{$key}{s_align} = "right";
		} elsif($$align{$key}{q_end} <= (1 + $dev) && $$align{$key}{s_start} <= (1 + $dev)){
		    $delete = 0;
		    $$align{$key}{reversed} = 1;
		    $$align{$key}{q_align} = "left";
		    $$align{$key}{s_align} = "left";
		}
	    }

	    #subject sequence reversed
	    if($$align{$key}{s_start} >= $$align{$key}{s_end}){
		if ($$align{$key}{q_end} >= ($$align{$key}{q_length} - $dev) && $$align{$key}{s_start} >= ($$align{$key}{s_length} - $dev)){
		    #both alignments at the end of sequences
		    $delete = 0;
		    $$align{$key}{reversed} = 2;
		    $$align{$key}{q_align} = "right";
		    $$align{$key}{s_align} = "right";
		}
		elsif($$align{$key}{q_start} <= (1 + $dev) && $$align{$key}{s_end} <= (1 + $dev)){
		    #both alignments on the right of sequences
		    $delete = 0;
		    $$align{$key}{reversed} = 2;
		    $$align{$key}{q_align} = "left";
		    $$align{$key}{s_align} = "left";
		}

	    }

	    #delete alignment when none of the cases is true
	    if($delete == 1){
		delete $align->{$key};
	    }
	}

	
	###filter double results (results for same alignment but switched query/subject)
	foreach my $key1 (keys %$align){
		foreach my $key2 (keys %$align){
		    unless (exists $$align{$key1} && exists $$align{$key2}){  #check if keys still exists
			next;
		    }

		    my $delete = 0;

		    if($$align{$key1}{subjectID} eq $$align{$key2}{queryID} 
		       && $$align{$key1}{queryID} eq $$align{$key2}{subjectID}){
			if ($$align{$key1}{reversed} == 0 && $$align{$key2}{reversed} == 0){
			    if ($$align{$key1}{s_start} == $$align{$key2}{q_start}
				&& $$align{$key1}{s_end} == $$align{$key2}{q_end}
				&& $$align{$key1}{q_start} == $$align{$key2}{s_start}
				&& $$align{$key1}{q_end} == $$align{$key2}{s_end}){
				$delete = 1;
			    }
			}else{
			    if ($$align{$key1}{s_start} == $$align{$key2}{q_end}
				&& $$align{$key1}{s_end} == $$align{$key2}{q_start}
				&& $$align{$key1}{q_start} == $$align{$key2}{s_end}
				&& $$align{$key1}{q_end} == $$align{$key2}{s_start}){
				$delete = 1;
			    }
			}
		    }

		    if($delete == 1){
			delete $align->{$key1};
		    }
		}
	}

	my $number = keys %$align;
	print LOG "First filtering step: search for overlap alignments\nresults: $number\n\n";
	return $align;
}


sub loadCoverage{
	my $bamfile = shift;

	###excluding bases with no coverage
  	#open (F, "samtools depth $bamfile | ") or die "Problem to open bam file $bamfile: $!\n";

	###including bases with no coverage
	open(F, "bedtools genomecov -d -ibam $bamfile | ") or die "Problem to open bam file $bamfile: $!\n";

 	my %coverage;

  	while (<F>){
		chomp;
		my ($chr,$pos,$depth)=split(/\t/);
		$coverage{$chr}{$pos}=$depth;
	}	
	return \%coverage;
}


sub evaluateOverlaps{
	my $overlaps = shift;
	my $coverage = shift;

	###calculate the mean coverage of the genome
	my $sum = 0;
	my $counter = 1;
	LOOP: foreach my $key1 (keys %$coverage){          #loop over sequences
	    foreach my $key2 (keys %{$$coverage{$key1}}){  #loop over position
		if ($counter > 1000000){                   #only use first million bases for calculation
		    last LOOP;
		}
		$sum = $sum + $$coverage{$key1}{$key2};
		$counter++;
	    }

	}
	my $mean = $sum / $counter;


	foreach my $key (keys %$overlaps){                 #loop over overlap alignments
		my $length = $$overlaps{$key}{alignment_length};
		my $s_length = $$overlaps{$key}{s_length};
		my $q_length = $$overlaps{$key}{q_length};
		my $s_id = $$overlaps{$key}{subjectID};
		my $q_id = $$overlaps{$key}{queryID};
		
		my $s_start;
		my $q_start;
		###get start base position of alignment
		if ($$overlaps{$key}{reversed} == 0){       #no reversed contig
		    $s_start = $$overlaps{$key}{s_start};
		    $q_start = $$overlaps{$key}{q_start};
		} elsif ($$overlaps{$key}{reversed} == 1){  #query reversed
		    $s_start = $$overlaps{$key}{s_start};
		    $q_start = $$overlaps{$key}{q_end}
		} elsif ($$overlaps{$key}{reversed} == 2){  #subject reversed
		    $s_start = $$overlaps{$key}{s_end};
		    $q_start = $$overlaps{$key}{q_start};
		}

		###calculate mean at overlap region
		my $sum_overlap = 0;
		for my $index (0..($length-1)){
			my $var1 = $s_start + $index;
			my $var2 = $q_start + $index;
			my $cov1 = $$coverage{$s_id}{$var1};
			my $cov2 = $$coverage{$q_id}{$var2}; 
	
			$sum_overlap = $sum_overlap + $cov1 + $cov2;
		}
		my $mean_overlap = $sum_overlap / ($length * 2);


		###testing if coverage at overlap is about half as high as overall coverage
		###if not: remove alignment
		my $frac = $mean_overlap / $mean;

		my $x = 0.1; #allowed deviation from 0.5
#		unless (0.5 - $x <= $frac && $frac <= 0.5 + $x){
		unless ($frac <= 0.5 + $x){
		    delete ($$overlaps{$key})
		}

	}

	return $overlaps;
}


sub getmerging{
    my $overlaps = shift;
    my $fasta = shift;


    ###get sequence names that are in the overlaps
    my %seqnames;
    ###get names of all querys/subjects in overlaps, stored as keys in sequences
    foreach my $key (keys %$overlaps){
	my $query = $$overlaps{$key}{queryID};
	my $subject = $$overlaps{$key}{subjectID};
	$seqnames{$query} = 1;
	$seqnames{$subject} = 1;
    }

    ### array to save the alignments that should be merged/not merged
    my @merging;
    my @nomerging;

    ###loop over contigs and check alignments
    for my $name (keys %seqnames){
	my %data;

	my @align;
	my @order;

	###loop over alignments for the current contig
	###check at which position the alignment is (left/right) and how long the merged sequence would be
	for my $key (keys %$overlaps){
	    if ($$overlaps{$key}{queryID} eq $name){
		push(@align, $key);
		push(@order, $$overlaps{$key}{q_align});
		# calculate length of the merged sequence
		my $length = $$overlaps{$key}{s_length} + $$overlaps{$key}{q_length} - 2 * $$overlaps{$key}{alignment_length};
		$data{$key} = {
		    length => $length,
		    position => $$overlaps{$key}{q_align}
		}
	    } elsif($$overlaps{$key}{subjectID} eq $name){
		push(@align, $key);
		push(@order, $$overlaps{$key}{s_align});
		# calculate length of the merged sequence
		my $length = $$overlaps{$key}{s_length} + $$overlaps{$key}{q_length} - 2 * $$overlaps{$key}{alignment_length};
		$data{$key} = {
		    length => $length,
		    position => $$overlaps{$key}{s_align}
		}
	    }
	}

	###determine how many alignments there are for the left and right end of the current contig
	my %counts;
	$counts{$_}++ for @order;

	###check if there are more than one alignment for one side of the sequence
	###decide which one to merge: choose alignment that results in longest merged sequence
	if ($counts{left} <= 1 && $counts{right} <= 1){
	    push(@merging, @align);    #if not more than one alignment for each side: all alignments can be merged
	    next;
	} else{
	    if ($counts{left} > 1 && $counts{right} <= 1){
		my $max = 0;
		my $maxkey;
		
		#get alignment on the left side that should be merged: the one with maximum length
		for my $i (keys %data){
		    if ($data{$i}{position} eq "left" && $data{$i}{length} > $max){
			$max = $data{$i}{length};
			$maxkey = $i;
		    }
		}
		
		push(@merging, $maxkey);

		#remember the ones that should not be merged
		for my $i (keys %data){
		    if ($data{$i}{position} eq "left" && $i ne $maxkey){
			push(@nomerging, $i);
		    }
		}	

		#save alignments on the right side in merging
		for my $i (keys %data){
		    if ($data{$i}{position} eq "right"){
			push(@merging, $i);
		    }
		}
	    } 
	    elsif ($counts{right} > 1 && $counts{left} <= 1){
		my $max = 0;
		my $maxkey;

		#get alignment on the right side that should be merged
		for my $i (keys %data){
		    if ($data{$i}{position} eq "right" && $data{$i}{length} > $max){
			$max = $data{$i}{length};
			$maxkey = $i;
		    }
		}
		push(@merging, $maxkey);

		#remember the ones that should not be merged
		for my $i (keys %data){
		    if ($data{$i}{position} eq "right" && $i ne $maxkey){
			push(@nomerging, $i);
		    }
		}
	       
		#save alignments on the left side in merging
		for my $i (keys %data){
		    if ($data{$i}{position} eq "left"){
			push(@merging, $i);
		    }
		}
	    }
	    else {
		###get alignment for left side
		my $max = 0;
		my $maxkey;
		
		#get alignment on the left side that should be merged: the one with maximum length
		for my $i (keys %data){
		    if ($data{$i}{position} eq "left" && $data{$i}{length} > $max){
			$max = $data{$i}{length};
			$maxkey = $i;
		    }
		}
		push(@merging, $maxkey);

		#remember the ones that should not be merged
		for my $i (keys %data){
		    if ($data{$i}{position} eq "left" && $i ne $maxkey){
			push(@nomerging, $i);
		    }
		}

		###get alignment for right side
		my $max = 0;
		my $maxkey;

		#get alignment on the right side that should be merged
		for my $i (keys %data){
		    if ($data{$i}{position} eq "right" && $data{$i}{length} > $max){
			$max = $data{$i}{length};
			$maxkey = $i;
		    }
		}
		push(@merging, $maxkey);

		#remember the ones that should not be merged
		for my $i (keys %data){
		    if ($data{$i}{position} eq "right" && $i ne $maxkey){
			push(@nomerging, $i);
		    }
		}
	    }	    
	}
    }

    ###remove duplicates from merging and nomerging
    my %hash1 = map { $_, 1} @merging;
    my @merging = keys %hash1;
    my %hash2 = map { $_, 1} @nomerging;
    my @nomerging = keys %hash2;

    ###get values that are in merging, but not in nomerging (stored in results)
    my %seen;     
    my @results;

    @seen{@nomerging} = ( );

    foreach my $item (@merging) {
	push(@results, $item) unless exists $seen{$item};
    }

    ###keep only alignments that should be merged
    
    my %results = map { $_ => 1 } @results;
    foreach my $key (keys %$overlaps){
	## if key is in results, then ok
	## if not: delete
	unless (exists($results{$key})){
	    delete $overlaps -> {$key};
	} 
    }

    my %newseqnames;
    ###get names of all querys/subjects in overlaps, stored as keys in sequences
    foreach my $key (keys %$overlaps){
	my $query = $$overlaps{$key}{queryID};
	my $subject = $$overlaps{$key}{subjectID};
	$newseqnames{$query} = 1;
	$newseqnames{$subject} = 1;
    }


    my $number = keys %$overlaps;
    print LOG "Third filtering step: get best alignments when there are more than one for one sequence \nalignments: $number\n";

    return $overlaps;
}


sub readseq{
    my $fastafile = shift;
    my $overlaps = shift;          #reference to hash
    my $contained = shift;         #reference to array
    my $notcovered = shift;        #reference to array
    my $newname = shift;
    my $containedname = shift;    
    my $notcoveredname = shift;

    my %sequences;
    ###get names of all querys/subjects in overlaps, stored as keys in sequences
    foreach my $key (keys %$overlaps){
	my $query = $$overlaps{$key}{queryID};
	my $subject = $$overlaps{$key}{subjectID};
	$sequences{$query} = 1;
	$sequences{$subject} = 1;
    }


    ###get names of all contained sequences, stored as keys in contained
    my %contained = map { $_ => 1 } @$contained;

    ###get names of all not covered sequences, stored as keys in notcovered
    my %notcovered = map { $_ => 1 } @$notcovered;

    ###extract sequences of all overlapping sequences from fasta file, stored in hash
    ###write non-overlap sequences in new fasta file
    ###write contained sequences in another new fasta file
    open(fasta, $fastafile) or die "problem opening fasta file $fastafile: $!\n";
    open(newfasta, ">$newname") or die "problem writing to new fasta file $newname: $!\n";
    open(containedfasta, ">$containedname") or die "problem writing to new fasta file $containedname: $!\n";
    open(notcoveredfasta, ">$notcoveredname") or die "problem writing to new fasta file $notcoveredname: $!\n";

    my %seqs;
    my $header;
    my $test = 0;
    
    while (<fasta>){
	chomp();
	if(/^>(.+)/){
	    my $name = $1;
	    if(exists($sequences{$name})){         #check for overlap sequences
		$test = 1;
		$header = $name;
		next;
	    } elsif(exists($contained{$name})){    #check for contained sequences
		$test = 2;
	    } elsif(exists($notcovered{$name})){   #check for not covered sequences
		$test = 3;
	    } else{
		$test = 4;
	    }
	}
	
	if($test == 1){
	    $seqs{$header} = $seqs{$header}.$_;    #save overlap sequences in hash
	} elsif ($test == 2){
	    print containedfasta "$_\n";           #write contained sequences in additional fasta file
	} elsif ($test == 3){
	    print notcoveredfasta "$_\n";          #write not covered sequences in additional fasta file
	} else {
	    print newfasta "$_\n";                 #write non- overlap sequences into new fasta file
	}
    }
    close(fasta);

    my $number = keys %seqs;
    print LOG "sequences: $number\n";

    return \%seqs;
}


sub mergeoverlaps{
    my $overlaps = shift;
    my $seq = shift;
    my $resultname = shift;

    my %mergedseq;

    open(results, ">$resultname") or die "problem writing result file $resultname: $!\n";

    ###get groups of sequences that merge together
    my $ref_groups = getgroups($overlaps);
    print Dumper $ref_groups;

    my $number = keys %$ref_groups;
    print LOG "merged sequences: $number\n\n";

    ###change structure of overlaps for better access to data
    my $hash = changestructure($overlaps);

    ###loop over mergings/groups
    foreach my $merge (keys %$ref_groups){
	my @names = keys %{$$ref_groups{$merge}};    #get all sequence names of this merge
	my $startseq;

	my $newname = "newseq".$merge;
	print results "$newname\n";


	###find a start sequence: sequence that is just in one alignment
	###if possible: sequence where alignment is on the right side

	STARTSEQ: foreach my $name (@names){
	    my %helper;
	    foreach my $key (keys %$hash){
		if (exists($$hash{$key}{$name})){
		    $helper{$name}{$key} = 1;
		}

	    }

	    if (keys %{$helper{$name}} == 1){
		$startseq = $name;
		my $key = (keys %{$helper{$name}})[0];
		if ($$hash{$key}{$name}{align} eq "right"){  #stop search when seq. with alignment on right side is found
		    last STARTSEQ;                           #otherwise continue and try to find another one
		}
	    }
	}

	if ($startseq eq ""){
	    print results "problem: no start sequence\n";
	    next;
	}

	###go through alignments in the order of the sequences

	my $firstseq = $startseq;
	my $secondseq;
	my $alignment;
	delete $$ref_groups{$merge}{$firstseq};

	my $mergedseq;

	###do merging as long as there are sequences left in this group
	while(keys %{$$ref_groups{$merge}} > 0){

	    foreach my $key (keys %$hash){                   #loop over alignments
		if(exists($$hash{$key}{$firstseq})){         #find alignment including the current sequence
		    $alignment = $key;                       #save alignment number
		    my @both = keys %{$$hash{$key}};         #get both sequence names of that alignment
		    for my $i (0,1){                         
			if ($both[$i] eq $firstseq){
			    splice(@both, $i, 1);            #remove the first sequence from that array
			}
		    }
		    $secondseq = @both[0];                   #get the second sequence
		}
	    }

print "$firstseq \t$$hash{$alignment}{$firstseq}{start} \t$$hash{$alignment}{$firstseq}{end}";
	    print results "$firstseq \t$$hash{$alignment}{$firstseq}{start} \t$$hash{$alignment}{$firstseq}{end} \t";
	    print results "$secondseq \t$$hash{$alignment}{$secondseq}{start} \t$$hash{$alignment}{$secondseq}{end} \n";

	    my $addseq;

	    my $reversed = $$overlaps{$alignment}{reversed};

	    ###START SEQUENCE
	    ###at the beginning: take start sequence completely (excuding bases right of the alignment) 
	    if ($firstseq eq $startseq){

		###case 1: alignment on the right for first sequence, no reversion in alignment
		if ($$hash{$alignment}{$firstseq}{align} eq "right" && $reversed == 0){
		    $mergedseq = substr($$seq{$firstseq}, 0, $$hash{$alignment}{$firstseq}{end});
		}
		###case 2: alignment on the right for first sequence, reversion
		###decide not to reverse first sequence, but reverse second sequence later
		elsif ($$hash{$alignment}{$firstseq}{align} eq "right" && $reversed != 0){
		    #check if first sequence was marked as reversed to determine the end position
		    my $end;
		    if($$hash{$alignment}{$firstseq}{reversed} == 0){
			$end = $$hash{$alignment}{$firstseq}{end};
		    }else{
			$end = $$hash{$alignment}{$firstseq}{start};
		    }
		    $mergedseq = substr($$seq{$firstseq}, 0, $end);
		}
		###case 3: alignment on the left for first sequence, no reversion 
		###(can only happen when there is no start seq with aligment on the right)
		###reverse this sequence, also reverse second one later
		elsif ($$hash{$alignment}{$firstseq}{align} eq "left" && $reversed == 0){
		    my $seq = substr($$seq{$firstseq}, $$hash{$alignment}{$firstseq}{start});
		    $mergedseq = rev_comp($seq);
		}
		###case 4: alignment on the left for first sequence, reversion
		###reverse this sequence, don't reverse the second one
		elsif ($$hash{$alignment}{$firstseq}{align} eq "left" && $reversed != 0){
		    #determine end position
		    my $end;
		    if($$hash{$alignment}{$firstseq}{reversed} == 0){
			$end = $$hash{$alignment}{$firstseq}{start};
		    }else{
			$end = $$hash{$alignment}{$firstseq}{end};
		    }
		    my $seq = substr($$seq{$firstseq}, $end);
		    $mergedseq = rev_comp($seq);
		}
	    }

	    
	    ###if not at the beginning: remove bases on the right end of the last alignment
	    if ($firstseq ne $startseq){

		###get number of bases
		my $number;
		###no reversion, alignment on the right
		if ($$hash{$alignment}{$firstseq}{align} eq "right" && $$hash{$alignment}{$firstseq}{reversed} == 0){
		    $number = $$hash{$alignment}{$firstseq}{length} - $$hash{$alignment}{$firstseq}{end};
		}
		###reversion, alignment on the right of first seq
		elsif ($$hash{$alignment}{$firstseq}{align} eq "right" && $$hash{$alignment}{$firstseq}{reversed} != 0){
		    $number = $$hash{$alignment}{$firstseq}{length} - $$hash{$alignment}{$firstseq}{start}
		}

		###no reversion, alignment on the left
		elsif ($$hash{$alignment}{$firstseq}{align} eq "left" && $$hash{$alignment}{$firstseq}{reversed} == 0){
		    $number = $$hash{$alignment}{$firstseq}{start};
		}
	    
		###reversion, alignment on the left of first seq
		elsif ($$hash{$alignment}{$firstseq}{align} eq "left" && $$hash{$alignment}{$firstseq}{reversed} != 0){
		    $number = $$hash{$alignment}{$firstseq}{end};
		}

		###remove these bases
		my $length = length($mergedseq);
		$mergedseq = substr($mergedseq, 0, $length - $number);
	    }


	    ###FOLLOWING SEQUENCE
	    ###add the second sequence without the overlap

	    ###no reversion, alignment on the left
	    ###simple case, take seq. as it is
	    if ($$hash{$alignment}{$secondseq}{align} eq "left" && $reversed == 0){
		$addseq = substr($$seq{$secondseq}, $$hash{$alignment}{$secondseq}{end});
	    }
	    ###reversion, alignment on the left of first seq
	    ###don't reverse second one, first one already reversed
	    elsif ($$hash{$alignment}{$secondseq}{align} eq "left" && $reversed != 0){
		my $pos;
		if($$hash{$alignment}{$secondseq}{reversed} == 0){
		    $pos = $$hash{$alignment}{$secondseq}{end};
		}else{
		    $pos = $$hash{$alignment}{$secondseq}{start};
		}
		$addseq = substr($$seq{$secondseq}, $pos);		    
	    }

	    ###no reversion, alignment on the right (left for first seq)
	    ###reverse this sequence
	    elsif ($$hash{$alignment}{$secondseq}{align} eq "right" && $reversed == 0){
		my $seq = substr($$seq{$secondseq}, 0, $$hash{$alignment}{$secondseq}{start});
		$addseq = rev_comp($seq);
	    }
	    
	    ###reversion, alignment on the right of first seq
	    ###reverse second sequence now
	    elsif ($$hash{$alignment}{$secondseq}{align} eq "right" && $reversed != 0){
		my $pos;
		if($$hash{$alignment}{$secondseq}{reversed} == 0){
		    $pos = $$hash{$alignment}{$secondseq}{start};
		}else{
		    $pos = $$hash{$alignment}{$secondseq}{end};
		}
		my $seq = substr($$seq{$secondseq}, 0, $pos);
		$addseq = rev_comp($seq);
	    }


	    $mergedseq = $mergedseq.$addseq;

	    delete $$hash{$alignment};                     #delete alignment that has been merged ######
	    delete $$ref_groups{$merge}{$secondseq};       #delete second sequence from group (to stop loop when all sequences are gone)
	    $firstseq = $secondseq;
	}

	$mergedseq{$newname} = $mergedseq;
    }
    
    return \%mergedseq;
}


sub changestructure{
    my $overlaps = shift;

    ###takes the overlap hash and changes it into a structure with better access to the data
    my %hash;

    foreach my $key (keys %$overlaps){
	$hash{$key} = {$$overlaps{$key}{queryID} => {start => $$overlaps{$key}{q_start},
						    end => $$overlaps{$key}{q_end},
		                                    align => $$overlaps{$key}{q_align},
		                                    a_length => $$overlaps{$key}{alignment_length},
		                                    length => $$overlaps{$key}{q_length}},
		       $$overlaps{$key}{subjectID} => {start => $$overlaps{$key}{s_start},
						      end => $$overlaps{$key}{s_end},
                                                      align => $$overlaps{$key}{s_align},
                                                      a_length => $$overlaps{$key}{alignment_length},
		                                      length => $$overlaps{$key}{s_length}}
	};
	
	if($$overlaps{$key}{reversed} == 1){
	    ###query reversed
	    $hash{$key}{$$overlaps{$key}{queryID}}{reversed} = 1;
	}elsif($$overlaps{$key}{reversed} == 2){
	    ###subject reversed
	    $hash{$key}{$$overlaps{$key}{subjectID}}{reversed} = 1;
	}
    }

    return \%hash;
}


sub getgroups{
    
    my $overlaps = shift;

    ###determine which sequences belong together, save them in groups in hash
    my $countalign = 0;
    my %groups;

    LOOP: foreach my $key (keys %$overlaps){
	my $name1 = $$overlaps{$key}{queryID};
	my $name2 = $$overlaps{$key}{subjectID};

	###check if sequences already assigned
	foreach my $keyalign (keys %groups){
	    if (exists($groups{$keyalign}{$name1}) || exists($groups{$keyalign}{$name2})){
		next LOOP;
	    }
	}

	$countalign++;

	###save names in hash
	$groups{$countalign} = {$name1 => 1,
		    $name2 => 1};

	my $check = 1;
	while($check != 0){
	    my $oldsize = keys %{$groups{$countalign}};
	    foreach my $key2 (keys %$overlaps){
		my $name3 = $$overlaps{$key2}{queryID};
		my $name4 = $$overlaps{$key2}{subjectID};	    
		if (exists($groups{$countalign}{$name3}) || exists($groups{$countalign}{$name4})){
		    $groups{$countalign}{$name3} = 1;
		    $groups{$countalign}{$name4} = 1;
		}
	    }
	    my $newsize = keys %{$groups{$countalign}};
	    $check = $oldsize - $newsize;
	}
    }

    return \%groups;
}


sub rev_comp {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}


sub writenewfasta{
    my $newseq = shift;
    my $newname = shift;

    open(newfasta, ">>$newname") or die "problem to write to new fasta file $newname: $!\n";

    ###save merged sequences to new fasta file
    foreach my $key (keys %$newseq){
	my $header = ">" . "$key";
	print newfasta "$header\n";

	my $sequence = $$newseq{$key};
	my $length = 70;

	while (my $chunk = substr($sequence, 0, $length, "")){
	    print newfasta "$chunk\n";
	}
    }
}
