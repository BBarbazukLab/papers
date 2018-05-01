#! /usr/bin/perl
# Victor Amin 2010

use strict;
use warnings;
no strict "refs";

use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my %options;
getopts('A:B:q:o:i:l:h', \%options);

if ($options{h} || !$options{q} || !$options{A} || !$options{B} || !$options{l}) {Getopt::Std->version_mess(); HELP_MESSAGE(\*STDERR)}
sub HELP_MESSAGE {
	my $fh = shift;
	print $fh "\nThis script compares two SAM alignments, outputs counts to STDOUT.\n";
	print $fh "\tOPTIONS:\n";
	print $fh "\t-l [read length] [required]\n";
	print $fh "\t-q [starting fastq] [required]\n";
	print $fh "\t-A [first SAM] [required]\n";
	print $fh "\t-B [second SAM] [required]\n";
	print $fh "\t-o [table output]\n";
	print $fh "\t-i [comma separated list of group index to output]\n";
	print $fh "\n";
	print $fh "Alignments are assembled into a multidimensional hash that keeps track of the number of alignments, whether the alignment is exact, and chromosome and position of the alignment. Only alignments in the top 'strata' are counted, meaning if there is 1 exact match and 3 inexact matches, the alignment count for that read will be 1. Strata is determined by subtracting the edit distance from the number of read nucleotides aligned. In the event of an exact match, this number is equal to the read length. The hash is duplicated for the second alignment. Once both alignments hashes are complete, reads are compared one at a time and counted depending on the relationship between the two alignments.\n";
	exit;
}

my $read_length = $options{l};
my $group_list = $options{i} ? $options{i} : 0;

my %group_outputs = ();
if ($group_list ne "0") {
	my @gl = split(",", $group_list);
	for my $l (@gl) {
		$group_outputs{$l} = "$l";
		open $group_outputs{$l}, ">${l}.txt" or die "\nCould not open ${l}.txt: $!\n";	
	}
}

open FASTQ, "<$options{q}" or die "\nTrouble opening FASTQ file: $!\n";
open SAM_A, "<$options{A}" or die "\nTrouble opening SAM file A: $!\n";
open SAM_B, "<$options{B}" or die "\nTrouble opening SAM file B: $!\n";
if ($options{o}) {open TABLE, ">$options{o}" or die "\nCannot write to table: $!\n"}

(my $sam_a_name = $options{A}) =~ s/.sam$//i;
(my $sam_b_name = $options{B}) =~ s/.sam$//i;

print STDERR "\nCounting FASTQ...\n";
my %reads;
my $reads_count;
# Fill hash with zeros
#my %empty = (
#	chrom => "",
#	pos => "",
#	count => 0,
#	exact => 0,
#	edits => "$read_length" 
#);

my %empty = (
	chrom => "",
	count => 0,
	exact => 0,
	edits => 9 
);

while (<FASTQ>) {
	if ($.%4 == 1) {
		if (!/^\@/) {die "\nMalformed FASTQ on line $..\n"}
		chomp;
		s/^\@//;
		$reads{$_}{A} = { %empty };
		$reads{$_}{B} = { %empty };
		$reads_count++;
	}
}
close FASTQ;

my $reads_count_check = scalar keys %reads;

if ($reads_count_check != $reads_count) {die "\nRead counts don't match -- read ids may not be unique.\n"}

print STDERR "\nCounting SAM A...\n";
while (<SAM_A>) {
	if (/^\@/) {next}
	chomp;
	my @cols = split;
	
	if (!exists $reads{$cols[0]}) {
		print STDERR "WARNING: Read does not exist in FASTQ.\n";
		next;
	}
	
	if ($cols[1] == 4) {
		next;
	} else {
		my $edit_distance = 0;
		for my $column (@cols) {
			if ($column =~ /^NM\:i\:(\d+)$/) {
				$edit_distance = $1;
			}
		}
		
		my $match_mismatch_count = 0;	
		my @mms = $cols[5] =~ m/(\d+)M/g;
		for my $mms (@mms) {
			$match_mismatch_count += $mms;
		}
		
		$edit_distance += ($read_length - $match_mismatch_count); # edit distance is mismatches, insertions, deletions, plus the difference between aligned length and total length
		
		if ($edit_distance == 0) {
			if ($reads{$cols[0]}{A}{exact} == 1) {
				$reads{$cols[0]}{A}{count}++;
			} else {
				$reads{$cols[0]}{A}{edits} = 0;
				$reads{$cols[0]}{A}{exact} = 1;
				$reads{$cols[0]}{A}{count} = 1;
				$cols[2] =~ /^(.+?)\|/;
				$reads{$cols[0]}{A}{chrom} = $1;
				#$reads{$cols[0]}{A}{pos} = $cols[3];
			}
		} else {
			if ($edit_distance < $reads{$cols[0]}{A}{edits}) {
				$reads{$cols[0]}{A}{edits} = $edit_distance;
				$reads{$cols[0]}{A}{count} = 1;
				$cols[2] =~ /^(.+?)\|/;
				$reads{$cols[0]}{A}{chrom} = $1;
				#$reads{$cols[0]}{A}{pos} = $cols[3];
			} elsif ($edit_distance == $reads{$cols[0]}{A}{edits}) {
				$reads{$cols[0]}{A}{count}++;
			}
		}
	}
}
close SAM_A; 

print STDERR "\nCounting SAM B...\n";
while (<SAM_B>) {
	if (/^\@/) {next}
	chomp;
	my @cols = split;
	
	if (!exists $reads{$cols[0]}) {
		print STDERR "WARNING: Read does not exist in FASTQ.\n";
		next;
	}
	
	if ($cols[1] == 4) {
		next;
	} else {
		my $edit_distance = 0;
		for my $column (@cols) {
			if ($column =~ /^NM\:i\:(\d+)$/) {
				$edit_distance = $1;
			}
		}
		
		my $match_mismatch_count = 0;	
		my @mms = $cols[5] =~ m/(\d+)M/g;
		for my $mms (@mms) {
			$match_mismatch_count += $mms;
		}
		
		$edit_distance += ($read_length - $match_mismatch_count); # edit distance is mismatches, insertions, deletions, plus the difference between aligned length and total length

		if ($edit_distance == 0) {
			if ($reads{$cols[0]}{B}{exact} == 1) {
				$reads{$cols[0]}{B}{count}++;
			} else {
				$reads{$cols[0]}{B}{edits} = 0;
				$reads{$cols[0]}{B}{exact} = 1;
				$reads{$cols[0]}{B}{count} = 1;
				$cols[2] =~ /^(.+?)\|/;
				$reads{$cols[0]}{B}{chrom} = $1;
				#$reads{$cols[0]}{B}{pos} = $cols[3];
			}
		} else {
			if ($edit_distance < $reads{$cols[0]}{B}{edits}) {
				$reads{$cols[0]}{B}{edits} = $edit_distance;
				$reads{$cols[0]}{B}{count} = 1;
				$cols[2] =~ /^(.+?)\|/;
				$reads{$cols[0]}{B}{chrom} = $1;
				#$reads{$cols[0]}{B}{pos} = $cols[3];
			} elsif ($edit_distance == $reads{$cols[0]}{B}{edits}) {
				$reads{$cols[0]}{B}{count}++;
			}
		}
	}
}
close SAM_B;

my $a_only_single_exact = 0;
my $a_only_single_inexact = 0;
my $a_only_multi_exact = 0;
my $a_only_multi_inexact = 0;
my $b_only_single_exact = 0;
my $b_only_single_inexact = 0;
my $b_only_multi_exact = 0;
my $b_only_multi_inexact = 0;
my $both_unaligned = 0;
my $both_single_exact_same = 0;
my $both_single_exact_diff = 0;
my $both_single_inexact_same = 0;
my $both_single_inexact_diff = 0;
my $both_inexact_diff_equal = 0;
my $both_inexact_diff_a_better = 0;
my $both_inexact_diff_b_better = 0;
my $both_multi_exact = 0;
my $both_multi_inexact = 0;
my $a_se_b_si = 0;
my $a_si_b_se = 0;
my $a_se_b_me = 0;
my $a_me_b_se = 0;
my $a_se_b_mi = 0;
my $a_mi_b_se = 0;
my $a_si_b_me = 0;
my $a_me_b_si = 0;
my $a_si_b_mi = 0;
my $a_mi_b_si = 0;
my $a_me_b_mi = 0;
my $a_mi_b_me = 0;

my %counts;

print STDERR "\nCounting cases...\n";
if ($options{o}) {print TABLE "ID,A_CHROM,A_POS,A_COUNT,A_EXACT_FLAG,B_CHROM,B_POS,B_COUNT,B_EXACT_FLAG\n"}		
# another sanity check:
$reads_count_check = 0;
for my $id (keys %reads) {
	my $a_chrom = $reads{$id}{A}{chrom};
	#my $a_pos = $reads{$id}{A}{pos};
	#my $a_chrom = 'A';
	my $a_pos = 0;
	my $a_count = $reads{$id}{A}{count};
	my $a_exact = $reads{$id}{A}{exact};
	my $a_edits = $reads{$id}{A}{edits};
	
	my $b_chrom = $reads{$id}{B}{chrom};
	#my $b_pos = $reads{$id}{B}{pos};
	#my $b_chrom = 'B';
	my $b_pos = 1;
	my $b_count = $reads{$id}{B}{count};
	my $b_exact = $reads{$id}{B}{exact};
	my $b_edits = $reads{$id}{B}{edits};
	
	#$a_chrom =~ /^(.+?)\|/;
	my $fusion_id = $a_chrom;
	
	if (!exists $counts{$fusion_id}) {
		$counts{$fusion_id}{both_unaligned} = 0;
		$counts{$fusion_id}{b_only_single_exact} = 0;
		$counts{$fusion_id}{b_only_single_inexact} = 0;
		$counts{$fusion_id}{a_only_single_exact} = 0;
		$counts{$fusion_id}{a_only_single_inexact} = 0;
		$counts{$fusion_id}{both_single_exact_same} = 0;
		$counts{$fusion_id}{both_single_exact_diff} = 0;
		$counts{$fusion_id}{both_single_inexact_same} = 0;
		$counts{$fusion_id}{both_single_inexact_diff} = 0;
		$counts{$fusion_id}{both_inexact_diff_equal} = 0;
		$counts{$fusion_id}{both_inexact_diff_b_better} = 0;
		$counts{$fusion_id}{both_inexact_diff_a_better} = 0;
		$counts{$fusion_id}{a_exact_b_inexact} = 0;
		$counts{$fusion_id}{b_exact_a_inexact} = 0;
	}
				
	if ($a_count == 0 && $b_count == 0) { # none
		$both_unaligned++;
		$counts{$fusion_id}{both_unaligned}++;
		if (exists $group_outputs{9}) {print {$group_outputs{9}} "$id\n"}
	} elsif ($a_count == 0) { # b only
		if ($b_count == 1 && $b_exact == 1) {
			$b_only_single_exact++;
			$counts{$fusion_id}{b_only_single_exact}++;
			if (exists $group_outputs{5}) {print {$group_outputs{5}} "$id\n"}
		} elsif ($b_count == 1 && $b_exact == 0) {
			$b_only_single_inexact++;
			$counts{$fusion_id}{b_only_single_inexact}++;
			if (exists $group_outputs{6}) {print {$group_outputs{6}} "$id\n"}
		} elsif ($b_count > 1 && $b_exact == 1) {
			$b_only_multi_exact++;
			if (exists $group_outputs{7}) {print {$group_outputs{7}} "$id\n"}
		} elsif ($b_count > 1 && $b_exact == 0) {
		 	$b_only_multi_inexact++;
		 	if (exists $group_outputs{8}) {print {$group_outputs{8}} "$id\n"}
		} else {die "\nImpossible case occurred.\n"}
	} elsif ($b_count == 0) { # a only
		if ($a_count == 1 && $a_exact == 1) {
			$a_only_single_exact++;
			$counts{$fusion_id}{a_only_single_exact}++;
			if (exists $group_outputs{1}) {print {$group_outputs{1}} "$id\n"}
		} elsif ($a_count == 1 && $a_exact == 0) {
			$a_only_single_inexact++;
			$counts{$fusion_id}{a_only_single_inexact}++;
			if (exists $group_outputs{2}) {print {$group_outputs{2}} "$id\n"}
		} elsif ($a_count > 1 && $a_exact == 1) {
			$a_only_multi_exact++;
			if (exists $group_outputs{3}) {print {$group_outputs{3}} "$id\n"}
		} elsif ($a_count > 1 && $a_exact == 0) {
		 	$a_only_multi_inexact++;
		 	if (exists $group_outputs{4}) {print {$group_outputs{4}} "$id\n"}
		} else {die "\nImpossible case occurred.\n"}
	} elsif ($a_count > 0 && $b_count > 0) { # both
		if ($a_count == 1 && $a_exact == 1 && $b_count == 1 && $b_exact == 1) {
			if ($a_pos ne '' && $b_pos ne '' && $a_pos == $b_pos && $a_chrom eq $b_chrom) {
				$both_single_exact_same++;
				$counts{$fusion_id}{both_single_exact_same}++;
				if (exists $group_outputs{10}) {print {$group_outputs{10}} "$id\n"}
			} else {
				$both_single_exact_diff++;
				$counts{$fusion_id}{both_single_exact_diff}++;
				if (exists $group_outputs{11}) {print {$group_outputs{11}} "$id\n"}
			}
		} elsif ($a_count == 1 && $a_exact == 0 && $b_count == 1 && $b_exact == 0) {
			if ($a_pos ne '' && $b_pos ne '' && $a_pos == $b_pos && $a_chrom eq $b_chrom) {
				$both_single_inexact_same++;
				$counts{$fusion_id}{both_single_inexact_same}++;
				if (exists $group_outputs{12}) {print {$group_outputs{12}} "$id\n"}
			} else {
				$both_single_inexact_diff++;
				$counts{$fusion_id}{both_single_inexact_diff}++;
				
				if ($a_edits == $b_edits) {
					$both_inexact_diff_equal++;
					$counts{$fusion_id}{both_inexact_diff_equal}++;
				} elsif ($a_edits > $b_edits) {
					$both_inexact_diff_b_better++;
					$counts{$fusion_id}{both_inexact_diff_b_better}++;
				} elsif ($b_edits > $a_edits) {
					$both_inexact_diff_a_better++;
					$counts{$fusion_id}{both_inexact_diff_a_better}++;
				} else {
					die "\nImpossible case occurred.\n";
				}				
				if (exists $group_outputs{13}) {print {$group_outputs{13}} "$id\n"}
			}
		} elsif ($a_count > 1 && $a_exact == 1 && $b_count > 1 && $b_exact == 1) {
			$both_multi_exact++;
			if (exists $group_outputs{14}) {print {$group_outputs{14}} "$id\n"}
		} elsif ($a_count > 1 && $a_exact == 0 && $b_count > 1 && $b_exact == 0) {
			$both_multi_inexact++;
			if (exists $group_outputs{15}) {print {$group_outputs{15}} "$id\n"}
		} elsif ($a_count == 1 && $a_exact == 1 && $b_count == 1 && $b_exact == 0) {
			$a_se_b_si++;
			$counts{$fusion_id}{a_exact_b_inexact}++;
			if (exists $group_outputs{16}) {print {$group_outputs{16}} "$id\n"}
		} elsif ($a_count == 1 && $a_exact == 0 && $b_count == 1 && $b_exact == 1) {
			$a_si_b_se++;
			$counts{$fusion_id}{b_exact_a_inexact}++;
			if (exists $group_outputs{17}) {print {$group_outputs{17}} "$id\n"}
		} elsif ($a_count == 1 && $a_exact == 1 && $b_count > 1 && $b_exact == 1) {
			$a_se_b_me++;
			if (exists $group_outputs{18}) {print {$group_outputs{18}} "$id\n"}
		} elsif ($a_count > 1 && $a_exact == 1 && $b_count == 1 && $b_exact == 1) {
			$a_me_b_se++;
			if (exists $group_outputs{19}) {print {$group_outputs{19}} "$id\n"}
		} elsif ($a_count == 1 && $a_exact == 1 && $b_count > 1 && $b_exact == 0) {
			$a_se_b_mi++;
			if (exists $group_outputs{20}) {print {$group_outputs{20}} "$id\n"}
		} elsif ($a_count > 1 && $a_exact == 0 && $b_count == 1 && $b_exact == 1) {
			$a_mi_b_se++;
			if (exists $group_outputs{21}) {print {$group_outputs{21}} "$id\n"}
		} elsif ($a_count == 1 && $a_exact == 0 && $b_count > 1 && $b_exact == 1) {
			$a_si_b_me++;
			if (exists $group_outputs{22}) {print {$group_outputs{22}} "$id\n"}
		} elsif ($a_count > 1 && $a_exact == 1 && $b_count == 1 && $b_exact == 0) {
			$a_me_b_si++;
			if (exists $group_outputs{23}) {print {$group_outputs{23}} "$id\n"}
		} elsif ($a_count == 1 && $a_exact == 0 && $b_count > 1 && $b_exact == 0) {
			$a_si_b_mi++;
			if (exists $group_outputs{24}) {print {$group_outputs{24}} "$id\n"}
		} elsif ($a_count > 1 && $a_exact == 0 && $b_count == 1 && $b_exact == 0) {
			$a_mi_b_si++;
			if (exists $group_outputs{25}) {print {$group_outputs{25}} "$id\n"}
		} elsif ($a_count > 1 && $a_exact == 1 && $b_count > 1 && $b_exact == 0) {
			$a_me_b_mi++;
			if (exists $group_outputs{26}) {print {$group_outputs{26}} "$id\n"}
		} elsif ($a_count > 1 && $a_exact == 0 && $b_count > 1 && $b_exact == 1) {
			$a_mi_b_me++;
			if (exists $group_outputs{27}) {print {$group_outputs{27}} "$id\n"}
		} else {die "\nImpossible case occurred.\n"}
	} else {die "\nImpossible case occurred.\n"}
	
	if ($options{o}) {
		if ($a_count > 1) {
			$a_chrom = '*';
			$a_pos = '*';
		}
		if ($b_count > 1) {
			$b_chrom = '*';
			$b_pos = '*';
		}
		print TABLE "$id,$a_chrom,$a_pos,$a_count,$a_exact,$b_chrom,$b_pos,$b_count,$b_exact\n";
	}
	$reads_count_check++;
}
	
if ($reads_count_check != $reads_count) {die "\nRead count has changed unexpectedly.\n"}

my $total = $a_only_single_exact + $a_only_single_inexact + $a_only_multi_exact + $a_only_multi_inexact + $b_only_single_exact + $b_only_single_inexact + $b_only_multi_exact + $b_only_multi_inexact + $both_unaligned + $both_single_exact_same + $both_single_exact_diff + $both_single_inexact_same + $both_single_inexact_diff + $both_multi_exact + $both_multi_inexact + $a_se_b_si + $a_si_b_se + $a_se_b_me + $a_me_b_se + $a_se_b_mi + $a_mi_b_se + $a_si_b_me + $a_me_b_si + $a_si_b_mi + $a_mi_b_si + $a_me_b_mi + $a_mi_b_me;

if ($total != $reads_count) {print STDERR "\nWARNING: Total counts did not add up.\n"}

#open BETTER, ">fusion-alignment-comparison-counts.csv" or die "Could not open inexact alignments comparison table: $!\n";

print "FUSION_ID,BOTH_EXACT,BOTH_INEXACT_EQUAL,BERLIN_ONLY_EXACT,C1674_ONLY_EXACT,BERLIN_EXACT_C1674_INEXACT,C1674_EXACT_BERLIN_INEXACT,BERLIN_ONLY_SINGLE_INEXACT,C1674_ONLY_SINGLE_INEXACT,BERLIN_INEXACT_BETTER,C1674_INEXACT_BETTER\n";
for my $fusion_id (keys %counts) {
	print "$fusion_id,$counts{$fusion_id}{both_single_exact_diff},$counts{$fusion_id}{both_inexact_diff_equal},$counts{$fusion_id}{a_only_single_exact},$counts{$fusion_id}{b_only_single_exact},$counts{$fusion_id}{a_exact_b_inexact},$counts{$fusion_id}{b_exact_a_inexact},$counts{$fusion_id}{a_only_single_inexact},$counts{$fusion_id}{b_only_single_inexact},$counts{$fusion_id}{both_inexact_diff_a_better},$counts{$fusion_id}{both_inexact_diff_b_better}\n";	
}


print STDERR "
\tREADS\t$reads_count
1\t${sam_a_name} ONLY SINGLE EXACT\t$a_only_single_exact
2\t${sam_a_name} ONLY SINGLE INEXACT\t$a_only_single_inexact
3\t${sam_a_name} ONLY MULTI EXACT\t$a_only_multi_exact
4\t${sam_a_name} ONLY MULTI INEXACT\t$a_only_multi_inexact
5\t${sam_b_name} ONLY SINGLE EXACT\t$b_only_single_exact
6\t${sam_b_name} ONLY SINGLE INEXACT\t$b_only_single_inexact
7\t${sam_b_name} ONLY MULTI EXACT\t$b_only_multi_exact
8\t${sam_b_name} ONLY MULTI INEXACT\t$b_only_multi_inexact
9\tBOTH UNALIGNED\t$both_unaligned
10\tBOTH SINGLE EXACT SAME POS\t$both_single_exact_same
11\tBOTH SINGLE EXACT DIFF POS\t$both_single_exact_diff
12\tBOTH SINGLE INEXACT SAME POS\t$both_single_inexact_same
13\tBOTH SINGLE INEXACT DIFF POS\t$both_single_inexact_diff
\t\t\tEQUAL\t$both_inexact_diff_equal
\t\t\t${sam_a_name} BETTER\t$both_inexact_diff_a_better
\t\t\t${sam_b_name} BETTER\t$both_inexact_diff_b_better
14\tBOTH MULTI EXACT\t$both_multi_exact
15\tBOTH MULTI INEXACT\t$both_multi_inexact
16\t${sam_a_name} SINGLE EXACT ${sam_b_name} SINGLE INEXACT\t$a_se_b_si
17\t${sam_a_name} SINGLE INEXACT ${sam_b_name} SINGLE EXACT\t$a_si_b_se
18\t${sam_a_name} SINGLE EXACT ${sam_b_name} MULTI EXACT\t$a_se_b_me
19\t${sam_a_name} MULTI EXACT ${sam_b_name} SINGLE EXACT\t$a_me_b_se
20\t${sam_a_name} SINGLE EXACT ${sam_b_name} MULTI INEXACT\t$a_se_b_mi
21\t${sam_a_name} MULTI INEXACT ${sam_b_name} SINGLE EXACT\t$a_mi_b_se
22\t${sam_a_name} SINGLE INEXACT ${sam_b_name} MULTI EXACT\t$a_si_b_me
23\t${sam_a_name} MULTI EXACT ${sam_b_name} SINGLE INEXACT\t$a_me_b_si
24\t${sam_a_name} SINGLE INEXACT ${sam_b_name} MULTI INEXACT\t$a_si_b_mi
25\t${sam_a_name} MULTI INEXACT ${sam_b_name} SINGLE INEXACT\t$a_mi_b_si
26\t${sam_a_name} MULTI EXACT ${sam_b_name} MULTI INEXACT\t$a_me_b_mi
27\t${sam_a_name} MULTI INEXACT ${sam_b_name} MULTI EXACT\t$a_mi_b_me
\tTOTAL\t$total\n
";

