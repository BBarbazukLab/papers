#! /usr/bin/perl
# Victor Amin 2010

use strict;
use warnings;

use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my %options;
getopts('uqh', \%options);

if ($options{h}) {Getopt::Std->version_mess(); HELP_MESSAGE(\*STDERR)}
sub HELP_MESSAGE {
	my $fh = shift;
	print $fh "\nThis script converts FASTQ files to FASTA. STDIN/STDOUT. Status to STDERR.\n";
	print $fh "\tOPTIONS:\n";
	print $fh "\t-u [print only unique sequences]\n";
	print $fh "\t-q [print quality score to defline]\n";
	print $fh "Defline for unique FASTA is of the form ID:COUNT. If both -u and -q are\nspecified, quality score with highest total quality is selected for output.\nQuality string is last element of defline, delimited by ' | '.\n\n";
	exit;
}

my $unique_flag = $options{u} ? $options{u} : 0;
my $quality_flag = $options{q} ? $options{q} : 0;

my $SEQ_MODE = 1;
my $QUAL_MODE = 2;
my $mode = 1;

my %uniq_counts; # sequence => count
my %uniq_qual_sum; # sequence => sum of quals
my %uniq_qual_best; # sequence => quals
my %sequences; # ident => sequence
my $lines = 0;
my $reads_count = 0;
my $ident;
my $sequence;
my $quality;

while (<>) {
    chomp;
    if (/^\@/ && $mode == $SEQ_MODE) {
		s/\@//;
		$ident = $_;
		$reads_count++;
    } elsif (/^\+/) {
		$mode = $QUAL_MODE;
    } elsif ($mode == $SEQ_MODE) {
		$sequence .= $_;
		$lines++; # in case of multiline FASTQ file
    } elsif ($mode == $QUAL_MODE) {
		$quality .= $_;
		$lines--;
		if ($lines == 0) { # end of a record
			$mode = $SEQ_MODE;
			if (!$unique_flag) {
				if ($quality_flag) {$ident .= " | $quality"}
				$sequences{$ident} = $sequence;
			} else {
				$uniq_counts{$sequence}++;
				if ($quality_flag) {
					my @quals = split(//, $quality);
					map ($_ = (ord($_)-64), @quals); # converts solexa 1.3+ ascii to integers
					my $summation = join('+', @quals);
					my $sum = eval($summation);
					if (!$uniq_qual_sum{$sequence} || $sum > $uniq_qual_sum{$sequence}) {
						$uniq_qual_sum{$sequence} = $sum;
						$uniq_qual_best{$sequence} = $quality;
					}
				}
			}
			$sequence = '';
			$quality = '';
		}
    } else {
		die "\nError reading file on line $..\n";
    }
}

my $i = 0;
if ($unique_flag) {
	for my $seq (keys %uniq_counts) {
		$i++;
		my $defline = "$i:$uniq_counts{$seq}";
		if ($quality_flag) {
			$defline .= " | $uniq_qual_best{$seq}";
		}
		print ">$defline\n$seq\n";	
	}
} else {
	for my $ident (keys %sequences) {
		$i++;
		print ">$ident\n$sequences{$ident}\n";
	}
}

print STDERR "\n$i FASTQ sequences converted to $reads_count FASTA records.\n";
