#! /usr/bin/perl
# Victor Amin 2010

use strict;
use warnings;

use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my %options;
getopts('f:p:s:h', \%options);

if ($options{h} || !$options{p} || !$options{f} || !$options{s})  {Getopt::Std->version_mess(); HELP_MESSAGE(\*STDERR)}
sub HELP_MESSAGE {
	my $fh = shift;
	print $fh "\nThis script calculates stats for coverage from a consensus pileup. Table to STDOUT.\n";
	print $fh "\tOPTIONS:\n";
	print $fh "\t-p <file> consensus pileup [required]\n";
	print $fh "\t-f <file> reference FASTA [required]\n";
	print $fh "\t-s <file> alignments SAM [required]\n";
	print $fh "\n\n";
	exit;
}

open FASTA, "<$options{f}" or die "Could not open FASTA: $!\n";
open PILEUP, "<$options{p}" or die "Could not open pileup:$!\n";
open SAM, "<$options{s}" or die "Could not open sam: $!\n";

# Added read of SAM file
# purpose: to store number of mapped reads and read length
# for RPKM calculation
# Martin McCrory, 12/6/2010
print STDERR "Counting mapped reads...\n";
my $mapped_reads = 0;
my $read_length = 0;
while (<SAM>) {
	if (/^@/) {next}
	my ($read_id, $flag, $rname, $pos, $mapq, $cigar, $mrnm, $mpos, $isize, $seq) = split("\t",$_);
	if ($flag != 4) {
		$mapped_reads++;
	}
	
	if ($. < 10000) { # autodetect read length
		my $length = length $seq;
		if ($read_length < $length) {$read_length = $length}
	}
}
close SAM;


my $header;
my %lengths;
while (<FASTA>) {
	chomp;
	if (/^>/) {
		s/^>//;
		$header = $_;
		$header =~ s/\|.+$//;
	} else {
		$lengths{$header} += length $_;
	}
}

print "REFERENCE_ID,RPKM,APN,SD,CV,COUNT\n"; # If you aligned against fusions, you REFERENCE_ID will be the fusion id. If you aligned against a genome reference, this field will be the chromosome id.
my $last_seqname;
my @coverage;
my $coverage;
while (<PILEUP>) { # uses consensus pileup
	my ($seqname, $pos, $ref, $depth, $qual, $qual2) = split;
	$seqname =~ s/\|.+$//;
	
	if (!$last_seqname) {
		$last_seqname = $seqname;
	}
	
	if ($last_seqname ne $seqname) { # restart count for each new chromosome
		$lengths{$last_seqname} = 0;
		$last_seqname =~ s/\|.+$//;
		
		my $mean = $coverage/$lengths{$seqname}; # apn 
		my $diff = $lengths{$seqname} - @coverage;
		my $count = int(($coverage/76) + 0.99);

		if ($diff > 0) {
			for my $i (1..$diff) {
				push @coverage, 0;
			}
		}
		
		my $squares = 0;
		for my $cov (@coverage) {
			$squares += (($cov - $mean) * ($cov - $mean));
		}
		
		my $sd = sqrt($squares/$lengths{$seqname});
		my $cv = $sd/$mean;
		
		# Modified by McCrory
		# 12/6/2010
		if ($read_length != 0) {
			my $reads_in_exon = $coverage/$read_length;
		}
		else {
			my $reads_in_exon = 0;
		}
		my $rpkm = (1000000000 * $mean) / ($mapped_reads); # (APN * 10^9)/Mapped Reads
		
		print "$last_seqname,$rpkm,$mean,$sd,$cv,$count\n";
		
		@coverage = ();
		$coverage = 0;
		$last_seqname = $seqname;
	}
	
	$coverage += $depth;
	push @coverage, $depth;
	
}
close PILEUP;

for my $seqname (keys %lengths) {
	if ($lengths{$seqname} != 0) {
		print "$seqname,0,0,0,0,0\n";
	}
}


