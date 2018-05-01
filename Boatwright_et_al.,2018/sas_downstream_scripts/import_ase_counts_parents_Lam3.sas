
libname trago '/ufrc/barbazuk/lboat/T.lamottei_T.crocifolius/Complete_data_set/align_reads_for_samcompare/sas_data';


%macro import_cnts (reads, rep, ref, alt) ;

proc import OUT=WORK.&reads&rep._2_&ref._&alt 
            DATAFILE="/ufrc/barbazuk/lboat/T.lamottei_T.crocifolius/Complete_data_set/align_reads_for_samcompare/ase_counts/ase_counts_&reads._&rep._2_&ref._&alt..csv"
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=20; 
run;

data ase_&reads&rep._&ref._&alt._L3;
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

%import_cnts (Croc, 6, Lam, Croc) ;
%import_cnts (Croc, 7, Lam, Croc) ;
%import_cnts (Croc, 8, Lam, Croc) ;

%import_cnts (Lam, 3, Lam, Croc) ;
%import_cnts (Lam, 4, Lam, Croc) ;
%import_cnts (Lam, 5, Lam, Croc) ;

data trago.ase_4_bayes_Croc_reads_LC_L3 ;
set ase_Croc6_Lam_Croc_L3 ase_Croc7_Lam_Croc_L3 ase_Croc8_Lam_Croc_L3 ;
run ;

data trago.ase_4_bayes_Lam_reads_LC_L3 ;
set ase_Lam3_Lam_Croc_L3 ase_Lam4_Lam_Croc_L3 ase_Lam5_Lam_Croc_L3 ;
run ;

