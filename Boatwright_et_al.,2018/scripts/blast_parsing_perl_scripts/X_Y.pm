package X_Y;
require Exporter;
@ISA = qw(Exporter);
@EXPORT_OK = qw(XY);#subroutines to export

{
	#CURRENTLY ONLY WORKING FOR WU_BLAST - WHERE query aligns either +/-, rather than NCBI where
	#subject aligns +/-
	#provide list of values in this order:
	#query_length, sub_length, query_begin, query_end, subject_begin, subject_end, dir, allowed_extra
	#returns 1 or more if 'bad x-y
	
	
	
	sub XY {
		#vars
		my ($query_length, $subject_length, $query_begin, $query_end, $subject_begin, $subject_end, $dir, $allowed_extra) = @_;
		#flags
		my $y_left = 0;
		my $y_right = 0;
		
		#bad y-right
		if ($dir eq "+"){
			$y_left = 1 if ( ($subject_begin >= $allowed_extra) && ($query_begin >= $allowed_extra) ); #floating left end cannot be longer than 50bp 
								   
			$y_right = 1 if ( ($subject_end <= ($subject_length - $allowed_extra)) && ($query_end <= ($query_length - $allowed_extra)) );
								   #floating right end of STC cannot be longer than 50bp 
								   #unless match is within#25 bp of BAC start
		
		}
		elsif ($dir eq "-"){
			$y_left = 1 if ( ($query_length - $query_begin >= $allowed_extra) && ($subject_begin >= $allowed_extra) );
			$y_right = 1 if ( ($subject_length - $subject_end >= $allowed_extra) &&  ($query_end >= $allowed_extra) );
		
		}
		#print "Y-left = $y_left\tY-right = $y_right\n";
		my $y_value = ($y_left + $y_right);
		return $y_value;
	}	





1  #return true
}

