#!/usr/bin/env perl

use strict;
use Getopt::Std;
# inputs:
# -p <putative frameshift regions, in bed format>
# -s <Sam alignment format of evidence (rnaseq)>
my %argopts;
if (! getopts('p:s:', \%argopts)) { die "invalid arguments"; }

# ends of read to ignore when finding frameshift
my $READ_END_IGNORE = 5;
my $FS_FACTOR_THR = 10; # num reads saying FS must be this number of times greater than those saying no FS
my $MERGE_NEARBY_FS = 1; # find FS's caused by assembly correction from Arrow/Pilon on phased diploid assembly, where reads from heterozygous haplotigs may be used to correct a single site leading to an indel


# read SAM file
# $sam{$scaffold}{$startpos}{$cigar}{COUNT} = count, {ENDPOS}, {CIGARR}
my %sam = ();
my $line;
open (SAMFILE, $argopts{'s'}) || die $!;
while (defined ($line=<SAMFILE>)) {
	chomp $line;
	(my $readid, my $flag, my $scaf, my $startpos, my $mapq, my $cigar, my $rest1, my $rest2, my $rest3, my $readseq, my $rest4, my $tags) = split(/\t/, $line, 12);
	if ($cigar eq "*") { next; }
	if (index($cigar,"S")!=-1 || index($cigar,"H")!=-1) { next; }
	if (index($cigar,"P")!=-1) { die "unexpected cigar P $line\n";  }
        if (index($cigar,"=")!=-1 || index($cigar,"X")!=-1) {   die "unexpected cigar =X $line\n";  }

	$sam{$scaf}{$startpos}{$cigar}{COUNT}++;
	my @cigar_arr = split(/([MIDN])/, $cigar);
	$sam{$scaf}{$startpos}{$cigar}{CIGARR} = \@cigar_arr;
	my $currpos = $startpos+0;
	for (my $cigar_idx=0; $cigar_idx<scalar @cigar_arr; $cigar_idx+=2) {
		if ($cigar_arr[$cigar_idx+1] eq "D" || $cigar_arr[$cigar_idx+1] eq "M" || $cigar_arr[$cigar_idx+1] eq "N") {
			$currpos += $cigar_arr[$cigar_idx]+0;
		}
		elsif ($cigar_arr[$cigar_idx+1] eq "I") {}
		else {
			print "Weird cigar 2: $scaf\t$startpos\t$cigar\t".(join("\t",@cigar_arr)).", cigar_idx=$cigar_idx $cigar_arr[$cigar_idx+1]\n";
		}
	}
	$sam{$scaf}{$startpos}{$cigar}{ENDPOS} = $currpos-1;
    	# for insertions, we need to know the insertion string, so store the read
	if (index($cigar,"I")!=-1) {  
		$sam{$scaf}{$startpos}{$cigar}{READS}{$readseq}++;
	}
}
print "Read Sam alignments for ".(scalar keys %sam)." scaffolds\n";
#my $scaf = "opera_scaffold_6068";
#print "Scaf $scaf, num alignments = ".(scalar keys %{$sam{$scaf}})."\n";
#foreach my $startpos (sort keys %{$sam{$scaf}}) {
#	print "opera_scaffold_6068\t$startpos\t".(join "\t",%{$sam{$scaf}{$startpos}})."\n";
#}
close SAMFILE;





# check if there are any frameshift in regions
# $putative_FS{$scaf}{$pos}{D|1}{NUM} = number of supporting reads, {EVID}{startpos|cigar} = 1
#       Frameshift codes: {D|1}, {D|2} = deletions of 1 or 2 bps
#                         {I|1|G} = insertion of 1 bp "G"
#                         {I|2|GC} = insertion of 2 bp "GC"
#                         {N} = no frameshift here
# Also fill in end position of each read
my %putative_FS = ();
foreach my $scaf (keys %sam) {
	foreach my $startpos (keys %{$sam{$scaf}}) {
		foreach my $cigar (keys %{$sam{$scaf}{$startpos}}) {
			my $cigar_arr = $sam{$scaf}{$startpos}{$cigar}{CIGARR};
			my $currpos = $startpos;
			my $read_currpos = 1;
			for (my $cigar_idx=0; $cigar_idx<scalar @$cigar_arr; $cigar_idx+=2) {
				if ($$cigar_arr[$cigar_idx+1] eq "I") {
					if ($startpos+$READ_END_IGNORE <= $currpos && $currpos <= $sam{$scaf}{$startpos}{$cigar}{ENDPOS}-$READ_END_IGNORE) {  # don't consider ends of read
						if (!defined $sam{$scaf}{$startpos}{$cigar}{READS} || scalar keys %{$sam{$scaf}{$startpos}{$cigar}{READS}}==0) { die; }
						foreach my $read_string (keys %{$sam{$scaf}{$startpos}{$cigar}{READS}}) {
							my $ins_string = substr($read_string, $read_currpos-1, $$cigar_arr[$cigar_idx]+0);
							if ($ins_string eq "") { die "Ins string empty: $scaf $startpos $cigar $read_string"; }
 							# Insertion evidences are specific to the insertion string
							$putative_FS{$scaf}{$currpos}{"I|$$cigar_arr[$cigar_idx]|$ins_string"}{NUM} += $sam{$scaf}{$startpos}{$cigar}{READS}{$read_string};
							$putative_FS{$scaf}{$currpos}{"I|$$cigar_arr[$cigar_idx]|$ins_string"}{EVID}{"$startpos|$cigar|$sam{$scaf}{$startpos}{$cigar}{COUNT}"} = 1;
						}
					}
					$read_currpos += $$cigar_arr[$cigar_idx]+0;
				}
				elsif ($$cigar_arr[$cigar_idx+1] eq "D") {
					if ($startpos+$READ_END_IGNORE <= $currpos && $currpos <= $sam{$scaf}{$startpos}{$cigar}{ENDPOS}-$READ_END_IGNORE) {  # don't consider ends of read
						$putative_FS{$scaf}{$currpos}{"D|$$cigar_arr[$cigar_idx]"}{NUM} += $sam{$scaf}{$startpos}{$cigar}{COUNT};
						$putative_FS{$scaf}{$currpos}{"D|$$cigar_arr[$cigar_idx]"}{EVID}{"$startpos|$cigar|$sam{$scaf}{$startpos}{$cigar}{COUNT}"} = 1;
					}
					$currpos += $$cigar_arr[$cigar_idx]+0;
                                }
				elsif ($$cigar_arr[$cigar_idx+1] eq "M") {
					$currpos += $$cigar_arr[$cigar_idx]+0;
					$read_currpos += $$cigar_arr[$cigar_idx]+0;
				}
				elsif ($$cigar_arr[$cigar_idx+1] eq "N" ) {
					$currpos += $$cigar_arr[$cigar_idx]+0;
				}
	

			}
		}
	}

}


# only keep frameshifts in specified regions
if (defined $argopts{'p'}) {
	# read putative frameshift regions
	# $frameshift_regions{$scaf}{$start|$end} = 1
	my %frameshift_regions = ();
	open (BEDFILE, $argopts{'p'}) || die $!;
	while (defined ($line=<BEDFILE>)) {
		chomp $line;
		(my $scaf, my $start, my $end) = split(/\t/, $line);
		$frameshift_regions{$scaf}{$start+0} = $end+0;
	}
	close BEDFILE;


	# go through frameshift list to delete those outside frameshift regions
	foreach my $scaf (sort keys %putative_FS) {
		my @fspos_sorted = (sort {$a <=> $b} keys %{$putative_FS{$scaf}});
		my $start_fs_idx = 0;
		my %keep_fspos = ();
		foreach my $regionstart (sort {$a <=> $b} keys %{$frameshift_regions{$scaf}}) {
			# advance $start_fs_idx till curr fs is after region start
			while ($start_fs_idx < scalar @fspos_sorted && $fspos_sorted[$start_fs_idx] < $regionstart) {
				$start_fs_idx++;
			}
			if ($start_fs_idx >= scalar @fspos_sorted) {
				last;
			}
			my $regionend = $frameshift_regions{$scaf}{$regionstart}+0;
			for (my $curr_fs_idx = $start_fs_idx; $curr_fs_idx < scalar @fspos_sorted; $curr_fs_idx++) {
				if ($fspos_sorted[$curr_fs_idx] > $regionend) { last; }
				if ($regionstart <= $fspos_sorted[$curr_fs_idx] && $fspos_sorted[$curr_fs_idx] <= $regionend) {
					$keep_fspos{$fspos_sorted[$curr_fs_idx]} = 1;
				}
			}
		}
		foreach my $fspos (keys %{$putative_FS{$scaf}}) {
			if (!defined $keep_fspos{$fspos}) {
				delete $putative_FS{$scaf}{$fspos};
			}
		}
	}
}



# go through frameshift list to count number of reads saying no frameshift
foreach my $scaf (sort keys %putative_FS) {
	my @fspos_sorted = (sort {$a <=> $b} keys %{$putative_FS{$scaf}});
	my $start_fs_idx = 0;

	foreach my $readstart (sort {$a <=> $b} keys %{$sam{$scaf}}) {
		foreach my $cigar (keys %{$sam{$scaf}{$readstart}}) {
			# advance $start_fs_idx till curr fs is after read
			while ($start_fs_idx < scalar @fspos_sorted && $fspos_sorted[$start_fs_idx] < $readstart) {
				$start_fs_idx++;
			}
			if ($start_fs_idx >= scalar @fspos_sorted) {
				last; 
			}
			# for this read, check start fs and subsequent ones 
			my $readend = $sam{$scaf}{$readstart}{$cigar}{ENDPOS};
			for (my $curr_fs_idx = $start_fs_idx; $curr_fs_idx < scalar @fspos_sorted; $curr_fs_idx++) {
				if ($fspos_sorted[$curr_fs_idx] > $readend) { last; }
				if ($fspos_sorted[$curr_fs_idx] < $readstart+$READ_END_IGNORE || $fspos_sorted[$curr_fs_idx] > $readend-$READ_END_IGNORE) { next; } # this cannot be outside this loop
				my $cigar_arr = $sam{$scaf}{$readstart}{$cigar}{CIGARR};
				my $curr_cigar_pos = $readstart;
				for (my $cigar_idx=0; $cigar_idx<scalar @$cigar_arr; $cigar_idx+=2) {
					if ($$cigar_arr[$cigar_idx+1] eq "M" || $$cigar_arr[$cigar_idx+1] eq "D" || $$cigar_arr[$cigar_idx+1] eq "N") {
						$curr_cigar_pos += $$cigar_arr[$cigar_idx]+0;
					}
					if ($curr_cigar_pos == $fspos_sorted[$curr_fs_idx]) {
						if ($$cigar_arr[$cigar_idx+1] eq "I") {
							last;
						}
					}
					elsif ($curr_cigar_pos > $fspos_sorted[$curr_fs_idx]) {
						if ($$cigar_arr[$cigar_idx+1] eq "M") {
							$putative_FS{$scaf}{$fspos_sorted[$curr_fs_idx]}{N}{NUM} += $sam{$scaf}{$readstart}{$cigar}{COUNT};
							$putative_FS{$scaf}{$fspos_sorted[$curr_fs_idx]}{N}{EVID}{"$readstart|$cigar|$sam{$scaf}{$readstart}{$cigar}{COUNT}"} = 1;
						}
						last;
					}
				}
			}
		} # end foreach cigar
		if ($start_fs_idx >= scalar @fspos_sorted) {
			last;
		}
	}
}



# correct errors due to assembly of heterozygous haplotigs, which result in 2 nearby insertions / deletions
# if: 1. deletions X1 and X2 are near (within 3 bp)
#     2. X1's deletion-supporting reads are not X2's deletion-supporting reads
#     3. |X1_deletion_supporting_reads| ~= |X2_deletion_supporting_reads| 
#     3. Let X1 be the deletion with more supporting reads. Let R = X2's deletion-supporting reads that are X1's no-FS reads. 
#        If |X1_deletion_supporting_reads| + |R| > thr . (|X1_no_FS_reads| - |R|
#        Then move R from X1_no_FS_reads into X1_deletion_supporting_reads, delete X2
if ($MERGE_NEARBY_FS==1) {

my $multi_fs_dist_thr = 3;
foreach my $scaf (sort keys %putative_FS) {
  my @fspos_array = sort {$a <=> $b} keys %{$putative_FS{$scaf}};
  my %fspos_to_delete = (); # store FS to be deleted (by merging with another FS)
  for (my $fs_idx=0; $fs_idx < scalar @fspos_array - 1; $fs_idx++) {
    # Get nearby frameshifts of the same type
    my $x1pos = $fspos_array[$fs_idx];
    my $x2pos = $fspos_array[$fs_idx+1];
    if (defined $fspos_to_delete{$x1pos} || defined $fspos_to_delete{$x2pos}) { next; }
    #print "x1pos = $x1pos, x2pos = $x2pos\n";
    if ($x2pos - $x1pos > $multi_fs_dist_thr) { next; }
    # get the best frameshift types for x1 and x2
    my $num_FS = 0;
    my $x1_best_FS;
    foreach my $fs (sort keys %{$putative_FS{$scaf}{$x1pos}}) {
      if ($fs eq "N") { next; }
      if ($putative_FS{$scaf}{$x1pos}{$fs}{NUM} > $num_FS) {
        $num_FS = $putative_FS{$scaf}{$x1pos}{$fs}{NUM};
        $x1_best_FS = $fs;
      }
    }
    $num_FS = 0;
    my $x2_best_FS;
    foreach my $fs (sort keys %{$putative_FS{$scaf}{$x2pos}}) {
      if ($fs eq "N") { next; }
      if ($putative_FS{$scaf}{$x2pos}{$fs}{NUM} > $num_FS) {
        $num_FS = $putative_FS{$scaf}{$x2pos}{$fs}{NUM};
        $x2_best_FS = $fs;
      }
    }
    #print "x1 best FS = $x1_best_FS, x2 best FS = $x2_best_FS\n";
    if (substr($x1_best_FS,0,3) ne substr($x2_best_FS,0,3)) { next; }
    #print "Same type!\n";

    # Let x1 be the FS with more reads
    #print "x1 num reads = $putative_FS{$scaf}{$x1pos}{$x1_best_FS}{NUM}, reads = ".(join(",",sort keys %{$putative_FS{$scaf}{$x1pos}{$x1_best_FS}{EVID}}))."\n";
    #print "x2 num reads = $putative_FS{$scaf}{$x2pos}{$x2_best_FS}{NUM}, reads = ".(join(",",sort keys %{$putative_FS{$scaf}{$x2pos}{$x2_best_FS}{EVID}}))."\n";
    if ($putative_FS{$scaf}{$x2pos}{$x2_best_FS}{NUM} > $putative_FS{$scaf}{$x1pos}{$x1_best_FS}{NUM}) {
      my $tmp_pos = $x2pos;
      my $tmp_best_FS = $x2_best_FS;
      $x2pos = $x1pos;
      $x2_best_FS = $x1_best_FS;
      $x1pos = $tmp_pos;
      $x1_best_FS = $tmp_best_FS;
      #print "SWITCHED, x1pos = $x1pos, x1_fs = $x1_best_FS, x1 num reads = $putative_FS{$scaf}{$x1pos}{$x1_best_FS}{NUM}, x2pos = $x2pos, x2_fs = $x2_best_FS, x2 num reads = $putative_FS{$scaf}{$x2pos}{$x2_best_FS}{NUM}\n"
    }
    # check that X1 and X2 have similar number of supporting reads
    if ($putative_FS{$scaf}{$x2pos}{$x2_best_FS}{NUM}/$putative_FS{$scaf}{$x1pos}{$x1_best_FS}{NUM}  < .5) {
      next;
    }
    #print "Number of reads are not too different!\n";

    # check intersection of X1 and X2's FS-supporting reads
    my $num_x1_x2_intersect = 0;
    foreach my $x1_read (keys %{$putative_FS{$scaf}{$x1pos}{$x1_best_FS}{EVID}}) {
      if (defined $putative_FS{$scaf}{$x2pos}{$x2_best_FS}{EVID}{$x1_read}) { $num_x1_x2_intersect++; }
    }
    #print "x1-x2 intesect evids = $num_x1_x2_intersect\n";

    # Find X2's FS-supporting-reads that are non-FS in X1. $RR{EVID}, $RR{NUM}
    if (!defined $putative_FS{$scaf}{$x1pos}{N}) { print "WARNING! X1 N not defined. This shouldn't happen!"; next; } # this shouldn't happen because x2 reads should be x1 N
    #print "x1 num non-FS reads = $putative_FS{$scaf}{$x1pos}{N}{NUM}, reads = ".(join(",",sort keys %{$putative_FS{$scaf}{$x1pos}{N}{EVID}}))."\n"; 
    my %RR = ();
    $RR{NUM} = 0;
    foreach my $x2evid (keys %{$putative_FS{$scaf}{$x2pos}{$x2_best_FS}{EVID}}) {
       if (defined $putative_FS{$scaf}{$x1pos}{N}{EVID}{$x2evid}) {
        $RR{EVID}{$x2evid} = 1;
        my @evidtoks = split(/\|/, $x2evid);
        $RR{NUM} += $evidtoks[2]+0;
      }
    }
    #print "RR NUM = $RR{NUM}. RR: ".join(',', sort keys %{$RR{EVID}})."\n";

    # If |X1_deletion_supporting_reads| + |R| > thr . (|X1_no_FS_reads| - |R|
    #print "If combine x1 and x2, x1 num reads = ".($putative_FS{$scaf}{$x1pos}{$x1_best_FS}{NUM} + $RR{NUM}).", x1 num non-FS reads = ".($putative_FS{$scaf}{$x1pos}{N}{NUM} - $RR{NUM})."\n";
    if ($putative_FS{$scaf}{$x1pos}{$x1_best_FS}{NUM} + $RR{NUM} >= 5 & $putative_FS{$scaf}{$x1pos}{$x1_best_FS}{NUM} + $RR{NUM} > $FS_FACTOR_THR * ($putative_FS{$scaf}{$x1pos}{N}{NUM} - $RR{NUM})) {
      #print "Combining x1fs and x2fs\n";
      foreach my $newevid (keys %{$RR{EVID}}) {
        $putative_FS{$scaf}{$x1pos}{$x1_best_FS}{EVID}{$newevid} = 1;
        delete $putative_FS{$scaf}{$x1pos}{N}{EVID}{$newevid};
      }
      $putative_FS{$scaf}{$x1pos}{$x1_best_FS}{NUM} += $RR{NUM}; 
      $putative_FS{$scaf}{$x1pos}{N}{NUM} -= $RR{NUM};
      delete $putative_FS{$scaf}{$x2pos};
      #print "new x1 num reads = $putative_FS{$scaf}{$x1pos}{$x1_best_FS}{NUM}, x1 num non-FS reads = $putative_FS{$scaf}{$x1pos}{N}{NUM}\n";
    }
  }
}
} # end MERGE_NEARBY_FS==1



print "Putative frameshifts found in ".(scalar keys %putative_FS)." scaffolds\n";
my $num_total_putative_FS = 0;
my $num_FS_accepted = 0;
my %accepted_FS = (); # $accepted_FS{$scaf|$fspos}{NUMYES}, {NUMNO}, {DIFF}
my @accepted_FS_diffs;
foreach my $scaf (sort keys %putative_FS) {
	foreach my $fspos (sort {$a <=> $b} keys %{$putative_FS{$scaf}}) {
		# either it has > 1 FS type (so there is definitely something other than N), or it has 1 FS type and it is not N. (N = no frameshift)
		if (scalar keys %{$putative_FS{$scaf}{$fspos}} > 1 || !defined $putative_FS{$scaf}{$fspos}{N}) {
			$num_total_putative_FS++;
		}
		my $num_noFS = 0;
		if (defined $putative_FS{$scaf}{$fspos}{N}) {
			$num_noFS += $putative_FS{$scaf}{$fspos}{N}{NUM};
		}
		my $num_FS = 0;
		my $best_FS;
		foreach my $fs (sort keys %{$putative_FS{$scaf}{$fspos}}) {
			if ($fs eq "N") { next; }
			if ($putative_FS{$scaf}{$fspos}{$fs}{NUM} > $num_FS) {
				$num_FS = $putative_FS{$scaf}{$fspos}{$fs}{NUM};
				$best_FS = $fs;
			}
		}
		if ($num_FS >= 5 && $num_FS > $FS_FACTOR_THR*$num_noFS) {
			$num_FS_accepted++;
			print "$scaf\t$fspos\t$best_FS\t$putative_FS{$scaf}{$fspos}{$best_FS}{NUM}\tEvids = ".(join(",",sort keys %{$putative_FS{$scaf}{$fspos}{$best_FS}{EVID}})).
				"\t$num_noFS\tEvids = ".(join(",",sort keys %{$putative_FS{$scaf}{$fspos}{N}{EVID}}))."\n";
#			$accepted_FS{"$scaf|$fspos"}{NUMYES} = $putative_FS{$scaf}{$fspos}{$best_FS}{NUM};
#			$accepted_FS{"$scaf|$fspos"}{NUMNO} = 0;
#			if (defined $putative_FS{$scaf}{$fspos}{N}) {
#				print "$scaf\t$fspos\tN\t$putative_FS{$scaf}{$fspos}{N}{NUM}\tEvids = ".(join(",",sort keys %{$putative_FS{$scaf}{$fspos}{N}{EVID}}))."\n";
#				$accepted_FS{"$scaf|$fspos"}{NUMNO} = $putative_FS{$scaf}{$fspos}{N}{NUM};
#			}
#			$accepted_FS{"$scaf|$fspos"}{DIFF} = $accepted_FS{"$scaf|$fspos"}{NUMYES} - $accepted_FS{"$scaf|$fspos"}{NUMNO};
#			push (@accepted_FS_diffs, $accepted_FS{"$scaf|$fspos"}{DIFF});
		}
	}
}
print "Num putative FS = $num_total_putative_FS, num FS accepted = $num_FS_accepted\n";

#@accepted_FS_diffs = sort @accepted_FS_diffs;
#print "Accepted FS, num reads support - num reads against, min = $accepted_FS_diffs[0], 1q = ".$accepted_FS_diffs[floor(.25*(scalar @accepted_FS_diffs))].
#							", median = ".$accepted_FS_diffs[floor(.5*(scalar @accepted_FS_diffs))].", 3q = ".$accepted_FS_diffs[floor(.75*(scalar @accepted_FS_diffs))].
#							", max = ".$accepted_FS_diffs[(scalar @accepted_FS_diffs) - 1]."\n";



#foreach my $scaf (sort keys %putative_FS) {
#	foreach my $fspos (sort {$a <=> $b} keys %{$putative_FS{$scaf}}) {
#		foreach my $fs (sort keys %{$putative_FS{$scaf}{$fspos}}) {
#			print "$scaf\t$fspos\t$fs\tNum reads = $putative_FS{$scaf}{$fspos}{$fs}{NUM}\tEvids = ".(join(",",sort keys %{$putative_FS{$scaf}{$fspos}{$fs}{EVID}}))."\n";
#		}
#	}
#}




