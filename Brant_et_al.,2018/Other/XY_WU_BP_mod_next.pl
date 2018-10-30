#!/sw/bin/perl5.12.3 -w
use strict;

BEGIN   {
        push (@INC, "/ufrc/barbazuk/lboat/scripts");
}

use BPlitenew;
use X_Y('XY');

my $usage = '
WU_BP.pl <BLAST REPORT>
';
die $usage unless (-e $ARGV[0]);

open (BLAST, "$ARGV[0]") or die "Can not open $ARGV[0]\n";

#for blastx, query alignment coordinates in NT, sub coords in aa,  Alignment length based on AA identities based on AAs
#take only top hit, hash HSPs, and process to get longest alignment


my $multiple_report = new BPlite::Multi(\*BLAST );
 #print "mult rep = $multiple_report\n";
 
 #for each BLAST report.....
 QUERY:	while(my $blast = $multiple_report->nextReport) {
     my $query = $blast->query;
	# print "#$query\n";
	 my ($qname, $qlength) = $query =~ /(\S+).+\((\d+,*\d+)\sletters/;
	# $qname =~ s/\,//g; #remove ',' after name
	# print "$qname $qlength\n";
	 
	 SUBJECT:while(my $subject = $blast->nextSbjct) {
         #print "query name = $qname\n";
		 
		 #next QUERY if $subject =~ /transpos|transcriptase|retro|gag|pol|gag-pol|polyprotein|mudr|rire|copia|poson|olyprotein|MITE|IS-element/i;
		 
		 #TIGR ATH specific
		 my ($sname) = $subject->name =~ />(\S+)/;
		 my $slength = $subject->length;
		 while(my $hsp = $subject->nextHSP) {
		 
		 	my $score =$hsp->score;
			my $bits = $hsp->bits;
			my $id = $hsp->percent;
			my $pval = $hsp->P;
			my $match = $hsp->match;
			my $pos = $hsp->positive;
			my $alig_len = $hsp->length;
			my $qb = $hsp->queryBegin;	  
			my $qe = $hsp->queryEnd; 	  
			my $sb = $hsp->sbjctBegin;	  
			my $se = $hsp->sbjctEnd; 	  
			my $q_alig = $hsp->queryAlignment; 
			my $sub_alig = $hsp->sbjctAlignment; 
			my $alig_string = $hsp->alignmentString;
			my $q_gaps = $hsp->queryGaps;	  
			my $s_gaps = $hsp->sbjctGaps;	  
			my $strand;
			my $alig_p;
			if ($qb <=$qe){
				$alig_p = $qe - $qb;
				$strand = "+";
			} elsif ($qe <=$qb){
				 $alig_p = $qb - $qe;
				 $strand = "-";
			}
			#my $allowed_overhang = (0.10 * $alig_len);
			my $cons = ($pos/$alig_len)*100;
		 	next if ($qname eq $sname);
			$alig_p =~ s/\,//g; $qlength =~ s/\,//g;
			my $alig_test = $alig_p/$qlength;
			# ($pval <= 1e-100) 
			#next unless (($id >= 95) && ($pval <= 1E-10) && ($alig_len >= 150));#((abs($se-$sb)/$slength >= 0.7) || (abs($qe-$qb)/$qlength >= 0.7))
			#my $xyFLAG = XY($qlength, $slength, $qb, $qe, $sb, $se, $strand, $allowed_overhang); 
			#next unless ($qlength >= 1000) && ($slength>=1000);# && ($xyFLAG == 0); #&& ($qlength <= 5000) && ($slength<=5000);
		 		print "$qname\t$qlength\t$sname\t$slength\t$score\t$pval\t$id\t$cons\t$alig_len\t$qb\t$qe\t$sb\t$se\t$strand\n" ;#if ( $id >= 65 && $alig_len >= 150);
	next QUERY;
		 }
		 
		 
     }
    
 }

#QUERYNAME, QUERYLEN, SUBNAME, SUBLEN, SCORE, PVAL, ID, CONS, ALIG_LEN, QB, QE, SB, SE

		
__END__

###################
#column definitions
###################

#################################
# parseblast column definitions #
#################################

 $hsp->score;
 $hsp->bits;
 $hsp->percent;
 $hsp->P;
 $hsp->match;
 $hsp->positive;
 $hsp->length;
 $hsp->queryBegin;      $hsp->qb;
 $hsp->queryEnd;        $hsp->qe;
 $hsp->sbjctBegin;      $hsp->sb;
 $hsp->sbjctEnd;        $hsp->se;
 $hsp->queryAlignment;  $hsp->qa;
 $hsp->sbjctAlignment;  $hsp->sa;
 $hsp->alignmentString; $hsp->as;
 $hsp->queryGaps;       $hsp->qg;
 $hsp->sbjctGaps;       $hsp->sg;









 1: query begin
 2: query end
 3: query name
 4: sbjct begin
 5: sbjct end
 6: sbjct name
 7: raw score
 8: bits (normalized score)
 9: E-value
10: P-value (really, also the E value)
11: percent identity
12: number of matches
13: number of positive scores (similarities)
14: length of alignment
15: length of query
16: length of sbjct
17: number of gaps in query alignment
18: number of gaps in sbjct alignment

END

