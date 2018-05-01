
libname trago '/ufrc/barbazuk/lboat/Old_World_New_BED/sas_data';

%macro import_cnts (reads, rep, ref, alt) ;

proc import OUT=WORK.&reads&rep._2_&ref._&alt 
            DATAFILE="/ufrc/barbazuk/lboat/Old_World_New_BED/trago_output/ase_counts/parents/ase_counts_&reads._&rep._2_&ref._&alt..csv"
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=20; 
run;

data ase_&reads&rep._&ref._&alt;
retain fusion_ID rep total_count ase_count ase &ref._count &alt._count ;
set &reads&rep._2_&ref._&alt;
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
rep=&rep;
hybrid="&reads";
rename fusion_ID = commonID ;
                drop BOTH_EXACT BOTH_INEXACT_EQUAL SAM_A_ONLY_EXACT SAM_B_ONLY_EXACT SAM_A_EXACT_SAM_B_INEXACT SAM_B_EXACT_SAM_A_INEXACT
                SAM_A_ONLY_SINGLE_INEXACT SAM_B_ONLY_SINGLE_INEXACT SAM_A_INEXACT_BETTER SAM_B_INEXACT_BETTER;

	run;

%mend ;

%import_cnts (Tdu, 1, Tdu, Tpo) ;
%import_cnts (Tdu, 2, Tdu, Tpo) ;
%import_cnts (Tdu, 3, Tdu, Tpo) ;

%import_cnts (Tdu, 1, Tdu, Tpr) ;
%import_cnts (Tdu, 2, Tdu, Tpr) ;
%import_cnts (Tdu, 3, Tdu, Tpr) ;

%import_cnts (Tpo, 1, Tdu, Tpo) ;
%import_cnts (Tpo, 2, Tdu, Tpo) ;
%import_cnts (Tpo, 3, Tdu, Tpo) ;

%import_cnts (Tpr, 1, Tdu, Tpr) ;
%import_cnts (Tpr, 2, Tdu, Tpr) ;
%import_cnts (Tpr, 3, Tdu, Tpr) ;

data trago.ase_4_bayes_Tdu_reads_Tdu_Tpo ;
set ase_Tdu1_Tdu_Tpo ase_Tdu2_Tdu_Tpo ase_Tdu3_Tdu_tpo;
run ;

data trago.ase_4_bayes_Tdu_reads_Tdu_Tpr ;
set ase_Tdu1_Tdu_Tpr ase_Tdu2_Tdu_Tpr ase_Tdu3_Tdu_tpr;
run ;

data trago.ase_4_bayes_Tpo_reads_Tdu_Tpo ;
set ase_Tpo1_Tdu_Tpo ase_Tpo2_Tdu_Tpo ase_Tpo3_Tdu_tpo;
run ;

data trago.ase_4_bayes_Tpr_reads_Tdu_Tpr ;
set ase_Tpr1_Tdu_Tpr ase_Tpr2_Tdu_Tpr ase_Tpr3_Tdu_tpr;
run ;


