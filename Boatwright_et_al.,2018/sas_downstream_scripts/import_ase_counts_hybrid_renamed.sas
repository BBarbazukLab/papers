
libname trago '/ufrc/barbazuk/lboat/T.lamottei_T.crocifolius/Complete_data_set/align_reads_for_samcompare/sas_data';


%macro import_cnts (reads, rep, ref, alt) ;

proc import OUT=WORK.&reads._&rep._2_&ref._&alt 
            DATAFILE="/ufrc/barbazuk/lboat/T.lamottei_T.crocifolius/Complete_data_set/align_reads_for_samcompare/ase_counts/ase_counts_&reads._&rep._2_&ref._&alt..csv" 
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

%import_cnts (Cast, 3141, Lam, Croc) ;
%import_cnts (Cast, 3134, Lam, Croc) ;
%import_cnts (Cast, 3113, Lam, Croc) ;
%import_cnts (Cast, 1364, Lam, Croc) ;
%import_cnts (Cast, 1335, Lam, Croc) ;
%import_cnts (Cast, 1311, Lam, Croc) ;
%import_cnts (Cast, 1053, Lam, Croc) ;
%import_cnts (Cast, 1044, Lam, Croc) ;
%import_cnts (Cast, 1014, Lam, Croc) ;
%import_cnts (Cast, 242, Lam, Croc) ;
%import_cnts (Cast, 234, Lam, Croc) ;
%import_cnts (Cast, 221, Lam, Croc) ;

data trago.ase_4_bayes_Cast2_reads_LC ;
set ase_Cast_242_Lam_Croc ase_Cast_234_Lam_Croc ase_Cast_221_Lam_Croc;
run ;

data trago.ase_4_bayes_Cast10_reads_LC ;
set ase_Cast_1053_Lam_Croc ase_Cast_1044_Lam_Croc ase_Cast_1014_Lam_Croc;
run ;

data trago.ase_4_bayes_Cast13_reads_LC ;
set ase_Cast_1335_Lam_Croc ase_Cast_1364_Lam_Croc ase_Cast_1311_Lam_Croc; 
run ;

data trago.ase_4_bayes_Cast31_reads_LC ;
set ase_Cast_3113_Lam_Croc ase_Cast_3134_Lam_Croc ase_Cast_3141_Lam_Croc; 
run ;
