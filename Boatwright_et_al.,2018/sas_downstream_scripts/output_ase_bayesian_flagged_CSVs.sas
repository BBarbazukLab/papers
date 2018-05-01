
libname trago "/ufrc/barbazuk/lboat/Old_World_New_BED/sas_data";

/* 
output ASE results
*/


/* All parents */
proc export data = trago.bayes_flag_sig_tdu_for_tdu_tpr
outfile='/ufrc/barbazuk/lboat/Old_World_New_BED/ase_results/bayes_flag_sig_TDU_for_UR.csv'
dbms=csv replace ;
run;

proc export data = trago.bayes_flag_sig_tpr_for_tdu_tpr
outfile='/ufrc/barbazuk/lboat/Old_World_New_BED/ase_results/bayes_flag_sig_TPR_for_UR.csv'
dbms=csv replace ;
run;


/* Tm parents */
proc export data = trago.bayes_flag_sig_tdu_for_tdu_tpo
outfile='/ufrc/barbazuk/lboat/Old_World_New_BED/ase_results/bayes_flag_sig_TDU_for_UO.csv'
dbms=csv replace ;
run;

proc export data = trago.bayes_flag_sig_tpo_for_tdu_tpo
outfile='/ufrc/barbazuk/lboat/Old_World_New_BED/ase_results/bayes_flag_sig_TPO_for_UO.csv'
dbms=csv replace ;
run;


/* Hybrids */
proc export data = trago.bayes_flag_sig_Tm_tdu_tpo
outfile='/ufrc/barbazuk/lboat/Old_World_New_BED/ase_results/bayes_flag_sig_Tm_for_UO.csv'
dbms=csv replace ;
run;

proc export data = trago.bayes_flag_sig_Tms_tdu_tpr
outfile='/ufrc/barbazuk/lboat/Old_World_New_BED/ase_results/bayes_flag_sig_Tms_for_UR.csv'
dbms=csv replace ;
run;

