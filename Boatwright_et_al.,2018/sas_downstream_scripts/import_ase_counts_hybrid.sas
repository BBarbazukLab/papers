
libname trago '/ufrc/barbazuk/lboat/Old_World_New_BED/sas_data';

%macro import_cnts (reads, rep, ref, alt) ;

proc import OUT=WORK.&reads._&rep._2_&ref._&alt 
            DATAFILE="/ufrc/barbazuk/lboat/Old_World_New_BED/trago_output/ase_counts/hybrids/ase_counts_&reads._&rep._2_&ref._&alt..csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=20; 
run;

proc print data=WORK.&reads._&rep._2_&ref._&alt (obs=10) ;
run;

data ase_&reads._&rep._&ref._&alt ;
retain fusion_ID rep total_count ase_count ase &ref._count &alt._count ;
set &reads._&rep._2_&ref._&alt ;
total_count= BOTH_EXACT + BOTH_INEXACT_EQUAL + SAM_A_EXACT_SAM_B_INEXACT +
	SAM_A_INEXACT_BETTER +SAM_A_ONLY_EXACT +SAM_A_ONLY_SINGLE_INEXACT +
	SAM_B_EXACT_SAM_A_INEXACT +SAM_B_INEXACT_BETTER +SAM_B_ONLY_EXACT +SAM_B_ONLY_SINGLE_INEXACT;
&ref._count=SAM_A_EXACT_SAM_B_INEXACT + SAM_A_INEXACT_BETTER +
	SAM_A_ONLY_EXACT +SAM_A_ONLY_SINGLE_INEXACT;
&alt._count=SAM_B_EXACT_SAM_A_INEXACT +SAM_B_INEXACT_BETTER +
	SAM_B_ONLY_EXACT + SAM_B_ONLY_SINGLE_INEXACT;
ase_count=&ref._count+&alt._count;
if ase_count>0 then ase=&ref._count/ase_count;
	else ase=.;
rep = &rep ;
hybrid = "&reads" ; 
rename fusion_ID = commonID ;
                drop BOTH_EXACT BOTH_INEXACT_EQUAL SAM_A_ONLY_EXACT SAM_B_ONLY_EXACT SAM_A_EXACT_SAM_B_INEXACT SAM_B_EXACT_SAM_A_INEXACT
                SAM_A_ONLY_SINGLE_INEXACT SAM_B_ONLY_SINGLE_INEXACT SAM_A_INEXACT_BETTER SAM_B_INEXACT_BETTER;

	run;

%mend ;

%import_cnts (Tm, 1, tdu, tpo) ;
%import_cnts (Tm, 2, tdu, tpo) ;
%import_cnts (Tm, 3, tdu, tpo) ;

%import_cnts (Tms, 1, tdu, tpr) ;
%import_cnts (Tms, 2, tdu, tpr) ;
%import_cnts (Tms, 3, tdu, tpr) ;

data trago.ase_4_bayes_Tm_reads_tdu_tpo ;
set ase_Tm_1_tdu_tpo ase_Tm_2_tdu_tpo ase_Tm_3_tdu_tpo;
run ;

data trago.ase_4_bayes_Tms_reads_tdu_tpr ;
set ase_Tms_1_tdu_tpr ase_Tms_2_tdu_tpr ase_Tms_3_tdu_tpr;
run ;

