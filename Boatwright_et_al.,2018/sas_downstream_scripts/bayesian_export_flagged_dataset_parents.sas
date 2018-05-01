libname trago '/ufrc/barbazuk/lboat/Old_World_New_BED/sas_data';

/* 
merge with flag 
export for bayesian
*/
/* proc contents data = trago.ase_bayesian_Tm_sbys_with_flag;run ; */

%macro add_flag (reads, alt) ;

proc sort data = trago.Bayes_flag_&reads._reads_tdu_&alt. ;
by commonID ;
proc sort data = trago.Ase_bayes_&reads._rds_tdu_&alt._sbys ;
by commonID ;
run;

data ase_bayes_&reads._tdu_&alt._flg;
merge trago.Ase_bayes_&reads._rds_tdu_&alt._sbys (in=in1) trago.Bayes_flag_&reads._reads_tdu_&alt. (in=in2);
by commonID ;
if in1;
run ;


%mend ;

%add_flag (tdu, tpo);
%add_flag (tpo, tpo);
%add_flag (tdu, tpr);
%add_flag (tpr, tpr);

data trago.ase_bayes_tdu_tdu_tpo_flag ;
set ase_bayes_tdu_tdu_tpo_flg;
rename
 tdu_count_rep1 = LINE_TOTAL_1
 tdu_count_rep2 = LINE_TOTAL_2
 tdu_count_rep3 = LINE_TOTAL_3
 tpo_count_rep1 = TESTER_TOTAL_1
 tpo_count_rep2 = TESTER_TOTAL_2
 tpo_count_rep3 = TESTER_TOTAL_3 ;

 drop  total_count_rep1
	   total_count_rep2
 	   total_count_rep3 ;
 run;



data trago.ase_bayes_tpo_tdu_tpo_flag ;
set ase_bayes_tpo_tdu_tpo_flg;
rename
 tdu_count_rep1 = LINE_TOTAL_1
 tdu_count_rep2 = LINE_TOTAL_2
 tdu_count_rep3 = LINE_TOTAL_3
 tpo_count_rep1 = TESTER_TOTAL_1
 tpo_count_rep2 = TESTER_TOTAL_2
 tpo_count_rep3 = TESTER_TOTAL_3 ;

 drop  total_count_rep1
	   total_count_rep2
 	   total_count_rep3 ;
 run;

data trago.ase_bayes_tdu_tdu_tpr_flag ;
set ase_bayes_tdu_tdu_tpr_flg;
rename
 tdu_count_rep1 = LINE_TOTAL_1
 tdu_count_rep2 = LINE_TOTAL_2
 tdu_count_rep3 = LINE_TOTAL_3
 tpr_count_rep1 = TESTER_TOTAL_1
 tpr_count_rep2 = TESTER_TOTAL_2
 tpr_count_rep3 = TESTER_TOTAL_3 ;

 drop  total_count_rep1
           total_count_rep2
           total_count_rep3 ;
 run;



data trago.ase_bayes_tpr_tdu_tpr_flag ;
set ase_bayes_tpr_tdu_tpr_flg;
rename
 tdu_count_rep1 = LINE_TOTAL_1
 tdu_count_rep2 = LINE_TOTAL_2
 tdu_count_rep3 = LINE_TOTAL_3
 tpr_count_rep1 = TESTER_TOTAL_1
 tpr_count_rep2 = TESTER_TOTAL_2
 tpr_count_rep3 = TESTER_TOTAL_3 ;

 drop  total_count_rep1
           total_count_rep2
           total_count_rep3 ;
 run;


proc export data =trago.ase_bayes_tdu_tdu_tpo_flag 
 outfile="/ufrc/barbazuk/lboat/Old_World_New_BED/empirical_bayesian_parents_input/ase_bayes_tdu_tdu_tpo_flag.csv"
DBMS = csv REPLACE;
run;

proc export data =trago.ase_bayes_tpo_tdu_tpo_flag 
 outfile="/ufrc/barbazuk/lboat/Old_World_New_BED/empirical_bayesian_parents_input/ase_bayes_tpo_tdu_tpo_flag.csv"
DBMS = csv REPLACE;
run;

proc export data =trago.ase_bayes_tdu_tdu_tpr_flag
 outfile="/ufrc/barbazuk/lboat/Old_World_New_BED/empirical_bayesian_parents_input/ase_bayes_tdu_tdu_tpr_flag.csv"
DBMS = csv REPLACE;
run;

proc export data =trago.ase_bayes_tpr_tdu_tpr_flag
 outfile="/ufrc/barbazuk/lboat/Old_World_New_BED/empirical_bayesian_parents_input/ase_bayes_tpr_tdu_tpr_flag.csv"
DBMS = csv REPLACE;
run;

