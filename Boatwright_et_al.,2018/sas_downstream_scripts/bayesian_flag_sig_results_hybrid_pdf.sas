
libname trago "/ufrc/barbazuk/lboat/Old_World_New_BED/sas_data";

ods listing ;
/* hybrids */
%macro flag (hybrid) ;

data trago.bayes_flag_sig_&hybrid;
    set trago.emp_&hybrid;
	if (q4_q025 ne " ") and (q4_q025 > 0.5 or q4_q975 < 0.5) then flag_q4_sig = 1; 
		else flag_q4_sig = 0;
	if (q5_q025 ne " ") and (q5_q025 > 0.5 or q5_q975 < 0.5) then flag_q5_sig = 1; 
		else flag_q5_sig = 0;
    if (q6_q025 ne " ") and (q6_q025 > 0.5 or q6_q975 < 0.5) then flag_q6_sig = 1; 
		else flag_q6_sig = 0;
	if flag_q4_sig = 1 and flag_q5_sig=1 and flag_q6_sig = 1 then flag_sig_&hybrid =1 ;
		else flag_sig_&hybrid =0;
	*keep commonID flag_: ;
	run ;


    run;
%mend ;

%flag (Tm_tdu_tpo) ;
%flag (Tms_tdu_tpr) ;

/* plot hybrids */
%macro plots (hybrid) ;

ods pdf file="/ufrc/barbazuk/lboat/Old_World_New_BED/reports/bayes_flag_sig_&hybrid..pdf";

title "Distribution of mean (*ESC*){unicode theta} for Emp_&hybrid";
proc sgplot data = trago.bayes_flag_sig_&hybrid ;
histogram q4_mean_theta / fillattrs=graphdata1 transparency=0.5;
density q4_mean_theta / lineattrs=graphdata1;
histogram q5_mean_theta / fillattrs=graphdata2 transparency=0.5;
density q5_mean_theta / lineattrs=graphdata2;
histogram q6_mean_theta / fillattrs=graphdata3 transparency=0.5;
density q6_mean_theta / lineattrs=graphdata3;
xaxis label="mean (*ESC*){unicode theta}";
run ;

ods pdf close ;

%mend ;

%plots (Tm_tdu_tpo) ;
%plots (Tms_tdu_tpr) ;

/* plot hybrids */
%macro plots (hybrid) ;


title "Distribution of mean (*ESC*){unicode theta} for Emp_&hybrid Sig all 3";
proc sgplot data =  trago.bayes_flag_sig_&hybrid ;
where flag_sig_&hybrid =1 ;
histogram q4_mean_theta / fillattrs=graphdata1 transparency=0.5;
density q4_mean_theta / lineattrs=graphdata1;
histogram q5_mean_theta / fillattrs=graphdata2 transparency=0.5;
density q5_mean_theta / lineattrs=graphdata2;
histogram q6_mean_theta / fillattrs=graphdata3 transparency=0.5;
density q6_mean_theta / lineattrs=graphdata3;
xaxis label="mean (*ESC*){unicode theta}";
run ;


%mend ;

%plots (Tm_tdu_tpo) ;
%plots (Tms_tdu_tpr) ;



