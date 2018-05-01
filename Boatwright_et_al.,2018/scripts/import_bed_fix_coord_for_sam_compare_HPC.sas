*libname trago 'S:\UFGI_Trago\sas_data';
libname trago '/ufrc/barbazuk/lboat/Old_World_New_BED/sas_data';

/*
bed files converted to csv in terminal window (using sed)
import trago bed files containing commonID
'fix' start and end so all starts are 0 

*/


%macro bed_in (ONE, TWO) ;

data WORK.bed_&ONE._&TWO    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile "/ufrc/barbazuk/lboat/Old_World_New_BED/references/&ONE.-&TWO._overlaps_WRT_orthologs.csv" delimiter = ','
 MISSOVER DSD lrecl=32767 ;
        informat consensedID $36. ;
        informat start best32. ;
        informat end best32. ;
        informat commonID $54. ;
        informat  b_VAR5 best32. ;
        informat strand $1. ;
        format consensedID $36. ;
        format start best32. ;
        format end best32. ;
        format commonID $54. ;
        format  b_VAR5 best32. ;
        format strand $1. ;
     input
                 consensedID $
                 start
                 end
                 commonID $
                 b_VAR5
                 strand $
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;


%mend ;

%bed_in (TDU, tpo) ;
%bed_in (TDU, tpr) ;
%bed_in (TPO, tdu) ;
%bed_in (TPR, tdu) ;


/* 'fix' the tdu coordinates in the bed file for sam compare between TDU and TPO */
proc contents data = bed_tdu_tpo ;
run;

data modify_bed_TDU_tpo ;
retain start_tdu_2 end_tdu_2  ;
set bed_tdu_tpo ;
if start > 0 then 
	do ;
	start_2 = (start - start) and end_2;
    end_2 = (end - start) ;
	end ;
	else if start = 0 then 
		do ;
		start_2 = start ;
		end_2 = end;
		end;
 run;

 data trago.TDU_tpo_bed_for_sam_compare ;
 retain commonID start_2 end_2 ;
 set modify_bed_TDU_tpo ;
rename start_2 = start  end_2= end ; 
keep commonID start_2 end_2 ;
run;

proc export data = trago.TDU_tpo_bed_for_sam_compare 
outfile ='/ufrc/barbazuk/lboat/Old_World_New_BED/references/TDU_tpo_bed_for_sam_compare.bed'
dbms=tab replace ;
run;

/* 'fix' the tdu coordinates in the bed file for sam compare between TDU and TPR */
proc contents data = bed_tdu_tpr ;
run;

data modify_bed_TDU_tpr ;
retain start_2 end_2  ;
set bed_tdu_tpr ;
if start > 0 then 
	do ;
	start_2 = (start - start) and end_2;
    end_2 = (end - start) ;
	end ;
	else if start = 0 then 
		do ;
		start_2 = start ;
		end_2 = end;
		end;
 run;

 data trago.TDU_tpr_bed_for_sam_compare ;
 retain commonID start_2 end_2 ;
 set modify_bed_TDU_tpr ;
rename start_2 = start  end_2= end ; 
keep commonID start_2 end_2 ;
run;

proc export data = trago.TDU_tpr_bed_for_sam_compare 
outfile ='/ufrc/barbazuk/lboat/Old_World_New_BED/references/TDU_tpr_bed_for_sam_compare.bed'
dbms=tab replace ;
run;


/* 'fix' the tdu coordinates in the bed file for sam compare between TPO and TDU */
proc contents data = bed_tpo_tdu ;
run;

data modify_bed_TPO_tdu ;
retain start_2 end_2  ;
set bed_tpo_tdu ;
if start > 0 then 
	do ;
	start_2 = (start - start) and end_2;
    end_2 = (end - start) ;
	end ;
	else if start = 0 then 
		do ;
		start_2 = start ;
		end_2 = end;
		end;
 run;

 data trago.TPO_tdu_bed_for_sam_compare ;
 retain commonID start_2 end_2 ;
 set modify_bed_TPO_tdu ;
rename start_2 = start  end_2= end ; 
keep commonID start_2 end_2 ;
run;

proc export data = trago.TPO_tdu_bed_for_sam_compare 
outfile ='/ufrc/barbazuk/lboat/Old_World_New_BED/references/TPO_tdu_bed_for_sam_compare.bed'
dbms=tab replace ;
run;


/* 'fix' the tdu coordinates in the bed file for sam compare between TPR and TDU */
proc contents data = bed_tpr_tdu ;
run;

data modify_bed_TPR_tdu ;
retain start_2 end_2  ;
set bed_tpr_tdu ;
if start > 0 then 
	do ;
	start_2 = (start - start) and end_2;
    end_2 = (end - start) ;
	end ;
	else if start = 0 then 
		do ;
		start_2 = start ;
		end_2 = end;
		end;
 run;

 data trago.TPR_tdu_bed_for_sam_compare ;
 retain commonID start_2 end_2 ;
 set modify_bed_TPR_tdu ;
rename start_2 = start  end_2= end ; 
keep commonID start_2 end_2 ;
run;

proc export data = trago.TPR_tdu_bed_for_sam_compare 
outfile ='/ufrc/barbazuk/lboat/Old_World_New_BED/references/TPR_tdu_bed_for_sam_compare.bed'
dbms=tab replace ;
run;



/* identify commonIDs that are very different in length ...... */
/*
data tdu ;
set A_tdu_tpo ;
rename  b_VAR5 = b_VAR5_tdu
 consensedID = consensedID_tdu
 end = end_tdu
 start = start_tdu
 strand = strand_tdu ;
 run ;

data tpo ;
set A_tpo_tdu ;
rename  b_VAR5 = b_VAR5_tpo
 consensedID = consensedID_tpo
 end = end_tpo
 start = start_tpo
 strand = strand_tpo ;
 run ;
 proc sort data = tdu ;
 by commonID ;
 proc sort data = tpo ;
 by commonID ;
 run;

 data tdu_plus_tpo ;
 merge tdu tpo ;
 by commonID ;
 run;

data tdu_plus_tpo_2 ;
retain commonID consensedID_tdu consensedID_tpo  start_tdu end_tdu start_tpo end_tpo strand_tdu strand_tpo ;
set tdu_plus_tpo ;
run;


data modify_bed_TDU_tpo ;
retain start_tdu_2 end_tdu_2  ;
set tdu_plus_tpo_2 ;
if start_tdu > 0 then 
	do ;
	start_tdu_2 = (start_tdu - start_tdu) and end_tdu_2;
    end_tdu_2 = (end_tdu - start_tdu) ;
	end ;
	else if start_tdu = 0 then 
		do ;
		start_tdu_2 = start_tdu ;
		end_tdu_2 = end_tdu;
		end;
 run;
data check ;
set modify_bed ;
diff_tdu  = (end_tdu - start_tdu) ;
diff_tpo  = (end_tpo - start_tpo) ;
too_much = abs(diff_tdu - diff_tpo) ;
if too_much ge 10 then flag_10 = 1;
else flag_10 = 0;
if too_much ge 500 then flag_500 = 1;
else flag_500 = 0;
if too_much ge 200 then flag_200 = 1;
else flag_200 = 0;
if too_much ge 500 then flag_500 = 1;
else flag_500 = 0;
if too_much ge 1 then flag_1 = 1;
else flag_1 = 0;
run;

proc freq data = check ;
tables flag_200 flag_500 flag_10 flag_1;
run;

proc export data  = check outfile='/ufrc/barbazuk/lboat/Old_World_New_BED/references/check_bed_files.csv'
dbms=csv replace ;
run;
*/



