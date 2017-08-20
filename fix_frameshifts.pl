#!/usr/bin/env perl

use strict;
use Getopt::Std;
# inputs:
# -i <frameshifts file (output from find_frameshifts.pl)>
# -s <scaffolds>
# -o <fixed scaffolds>
my %argopts;
if (! getopts('i:s:o:', \%argopts)) { die "invalid arguments"; }



# read frameshifts file
my %frameshifts = ();
my $line;
open (FSFILE, $argopts{'i'}) || die $!;
while (defined ($line=<FSFILE>)) {
	chomp $line;
	(my $scaf, my $fspos, my $fs, my $numreads, my $rest) = split (/\t/, $line, 5);
	if (substr($fs,0,1) eq "D" || substr($fs,0,1) eq "I") {
		$fspos = $fspos-1; # change to base 0
		if (defined $frameshifts{$scaf}{$fspos}) { die "$scaf $fspos already exists!"; }
		$frameshifts{$scaf}{$fspos} = $fs;
	}
}
print "Frameshifts read for ".(scalar keys %frameshifts)." scaffolds\n";


# read assembly
open (SCAFFILE, $argopts{'s'}) || die $!;
my %scaffolds = ();
my $curr_scafid = "";
my @curr_scaf_seqlines = ();
while (defined ($line=<SCAFFILE>)) {
	chomp $line;
	if (substr($line,0,1) eq ">") {
		# save the previous seq
		if ($curr_scafid ne "") {
			$scaffolds{$curr_scafid}{S} = join("",@curr_scaf_seqlines);
		}

		($curr_scafid, my $rest) = split(/\s/, $line, 2);
		$curr_scafid = substr($curr_scafid,1);
		@curr_scaf_seqlines = ();
		$scaffolds{$curr_scafid}{H} = $line;
		$scaffolds{$curr_scafid}{S} = "";
		$scaffolds{$curr_scafid}{N} = substr($curr_scafid,15)+0;
	}
	else {
		if ($curr_scafid eq "") { die; }
		push (@curr_scaf_seqlines, $line);
	}
}
if ($curr_scafid ne "") {
	$scaffolds{$curr_scafid}{S} = join("",@curr_scaf_seqlines);
}
print "Num scaffolds read = ".(scalar keys %scaffolds)."\n";



# fix scaffolds
foreach my $scaf (keys %frameshifts) {
	print "$scaf, orig len = ".(length($scaffolds{$scaf}{S}))."\n";
	my $scaf_edit_offset = 0;
	foreach my $fspos (sort {$a <=> $b} keys %{$frameshifts{$scaf}}) {
		print "$scaf, curr len = ".(length($scaffolds{$scaf}{S})).", curr edit offset = $scaf_edit_offset\t";
		my $fs = $frameshifts{$scaf}{$fspos};
		(my $fstype, my $fslen, my $ins_string) = split(/\|/, $fs);
		print "FS: $fspos $fs $fstype $fslen $ins_string";
		if ($fstype eq "I") {
			if (length($ins_string) != $fslen) { die "$scaf $fspos $fs"; }
			#print " ins_string=$ins_string";
			substr($scaffolds{$scaf}{S}, $fspos+$scaf_edit_offset, 0, $ins_string);
			$scaf_edit_offset += $fslen;
		}
		elsif ($fstype eq "D") {
			substr($scaffolds{$scaf}{S}, $fspos+$scaf_edit_offset, $fslen, "");
			$scaf_edit_offset -= $fslen;
		}
		else { die "weird FS code $fs"; }
		print "\n";
		
	}
}


# print scaffolds
open (OUTPUTFILE, ">$argopts{'o'}") || die $!;
foreach my $scaf (sort {$scaffolds{$a}{N} <=> $scaffolds{$b}{N}} keys %scaffolds) {
	print OUTPUTFILE "$scaffolds{$scaf}{H}\n$scaffolds{$scaf}{S}\n";

}



