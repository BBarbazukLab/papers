#!/usr/bin/perl -w

use strict;
use warnings;

open FH, "<grouped_hits" or die "Couldn't open file";

my @contigs=();
my @temp = ();

while(<FH>){
	@temp = split(/\s+/,$_);
	push(@contigs,@temp);
}# end while

close FH;

my $fasta = $ARGV[0];

open NEW, "<$fasta" or die "Couldn't open file";

my @fasta_contigs=();

while(<NEW>){
	if($_ =~ m/>(\w+_Contig_\d+)\s/){
		push(@fasta_contigs,$1);
	}# end if
}# end while

close NEW;
open OUT, ">non-unique_contigs";

for my $seq (@fasta_contigs){
	if($seq ~~ @contigs){
		print OUT ">$seq\n";	
	}else{
		print ">$seq\n";	
	}# end if-else
}# end for

