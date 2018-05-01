#! /usr/bin/perl

use strict;
use warnings;

use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my %options;
getopts('b:s:h', \%options);

if ($options{h} || !$options{b} || !$options{s}) {Getopt::Std->version_mess(); HELP_MESSAGE(\*STDERR)}
sub HELP_MESSAGE {
	my $fh = shift;
	print $fh "\nThis program filters out SAM alignments that don't fall within features in a BED.\n";
	print $fh "\tOPTIONS:\n";
	print $fh "\t-b <file> BED features to filter by [required]\n";
	print $fh "\t-s <file> SAM alignments to filter [required]\n";
	print $fh "\n";
	exit;
}

open SAM, "$options{s}" or die "Could not open SAM file: $!\n";
open BED, "$options{b}" or die "Could not open BED file: $!\n";


print STDERR "Loading BED...\n";
my %coverage;
while (<BED>) {
	if (/^track\s/) {next}
	chomp;
	my ($chrom,$chrom_start,$chrom_end,$ident,$score,$strand,$t_start,$t_end,$rgb,$block_count,$block_sizes,$block_starts) = split;
	$chrom_start++; # 0-based to 1-based
	my @sizes = $block_sizes ? split(',', $block_sizes) : ();
	my @starts = $block_starts ? split(',', $block_starts) : ();
	if (@sizes > 1) {
		for my $i (0..$#sizes) {
			my $start = $chrom_start + $starts[$i];
			my $end = $chrom_start + $starts[$i] + $sizes[$i] - 1;
			
			for my $j ($start..$end) {
				$coverage{$chrom}{$j}++;
			}
		}
	} else {
		for my $j ($chrom_start..$chrom_end) {
			$coverage{$chrom}{$j}++;
		}
	}
}

print STDERR "Filtering SAM...\n";
while (<SAM>) {
	if (/^@/) {next}
	my ($qname, $flag, $rname, $pos, $mapq, $cigar, $mrnm, $mpos, $isize, $seq, $qual) = split('\t');
	my $length = length $qual;
	my $end = $pos + $length - 1;
	
	for my $i ($pos..$end) {
		if (exists $coverage{$rname}{$i}) {
			print; last;
		}	
	}
}
