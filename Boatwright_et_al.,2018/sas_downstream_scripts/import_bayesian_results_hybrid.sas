/********************************************************************************
* This script imports the output from the Bayesian Machine
********************************************************************************/

libname trago "/ufrc/barbazuk/lboat/Old_World_New_BED/sas_data";

%macro import_bayesian (hybrid) ;

     data WORK.Emp_&hybrid   ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile "/ufrc/barbazuk/lboat/Old_World_New_BED/empirical_bayesian_hybrids_output/PG_emp_bayesian_results_&hybrid._hybrid.csv"
	delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat commonID $54. ;
        informat q4_mean_theta best32. ;
        informat q4_q025 best32. ;
        informat q4_q975 best32. ;
        informat q5_mean_theta best32. ;
        informat q5_q025 best32. ;
        informat q5_q975 best32. ;
        informat q6_mean_theta best32. ;
        informat q6_q025 best32. ;
        informat q6_q975 best32. ;
        format commonID $54. ;
        format q4_mean_theta best12. ;
        format q4_q025 best12. ;
        format q4_q975 best12. ;
        format q5_mean_theta best12. ;
        format q5_q025 best12. ;
        format q5_q975 best12. ;
        format q6_mean_theta best12. ;
        format q6_q025 best12. ;
        format q6_q975 best12. ;
     input
                 commonID $
                 q4_mean_theta
                 q4_q025
                 q4_q975
                 q5_mean_theta
                 q5_q025
                 q5_q975
                 q6_mean_theta
                 q6_q025
                 q6_q975
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;
%mend ;

/*
%macro import_bayesian (reads, ref1, ref2) ;

     data WORK.Emp_&reads._for_&ref1._&ref2   ;
     %let _EFIERR_ = 0; 
     infile "/ufrc/barbazuk/lboat/T.lamottei_T.crocifolius/Complete_data_set/align_reads_for_samcompare/empirical_bayesian_parents_output/PG_emp_bayesian_results_&reads._parents.csv"
	delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat commonID $54. ;
        informat q4_mean_theta best32. ;
        informat q4_q025 best32. ;
        informat q4_q975 best32. ;
        informat q5_mean_theta best32. ;
        informat q5_q025 best32. ;
        informat q5_q975 best32. ;
        informat q6_mean_theta best32. ;
        informat q6_q025 best32. ;
        informat q6_q975 best32. ;
        format commonID $54. ;
        format q4_mean_theta best12. ;
        format q4_q025 best12. ;
        format q4_q975 best12. ;
        format q5_mean_theta best12. ;
        format q5_q025 best12. ;
        format q5_q975 best12. ;
        format q6_mean_theta best12. ;
        format q6_q025 best12. ;
        format q6_q975 best12. ;
     input
                 commonID $
                 q4_mean_theta
                 q4_q025
                 q4_q975
                 q5_mean_theta
                 q5_q025
                 q5_q975
                 q6_mean_theta
                 q6_q025
                 q6_q975
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  
     output WORK.Emp_&reads._for_&ref1._&ref2;
     run;
%mend ;
*/

%import_bayesian (Tm_tdu_tpo) ;
%import_bayesian (Tms_tdu_tpr) ;

%macro setdata (hybrid) ;

data trago.Emp_&hybrid;
     set Emp_&hybrid;
run;

%mend ;

%setdata (Tm_tdu_tpo) ;
%setdata (Tms_tdu_tpr) ;
