
/* scatter plot of q5_mean theta from TM vs TDU */
proc sort data = emp_tm;
by commonID ;
proc sort data = Emp_tdu_for_tdu_tpo;
by commonID ;
run ;

data tm ;
set emp_tm ;
keep commonID q5_mean_theta ;
rename q5_mean_theta = q5_mean_theta_Tm;
run;
data tdu ;
set Emp_tdu_for_tdu_tpo ;
keep commonID q5_mean_theta ;
rename q5_mean_theta = q5_mean_theta_Tdu;
run;


data Tm_Tdu_compare_q5  missing;
merge Tm (in=in1) tdu (in=in2) ;
by commonid ;
if in1 and in2 then output Tm_Tdu_compare_q5 ;
else output missing ;
run;  /* 117 in missing..... */


title 'plot of Tm q5_mean_theta vs Tdu q5_mean_theta';
proc sgplot data = Tm_Tdu_compare_q5 ;
scatter x = q5_mean_theta_Tdu y = q5_mean_theta_Tm ;
run ;


/* scatter plot of q5_mean theta from TMS vs TDU */
proc sort data = emp_tms;
by commonID ;
proc sort data = Emp_tdu_for_tdu_tpr;
by commonID ;
run ;

data tms ;
set emp_tms ;
keep commonID q5_mean_theta ;
rename q5_mean_theta = q5_mean_theta_Tms;
run;
data tdu_tpr ;
set Emp_tdu_for_tdu_tpr ;
keep commonID q5_mean_theta ;
rename q5_mean_theta = q5_mean_theta_Tdu;
run;


data Tms_Tdu_compare_q5  missing;
merge Tms (in=in1) tdu_tpr (in=in2) ;
by commonid ;
if in1 and in2 then output Tms_Tdu_compare_q5 ;
else output missing ;
run;  /* 117 in missing..... */


title 'plot of Tms q5_mean_theta vs Tdu q5_mean_theta';
proc sgplot data = Tms_Tdu_compare_q5 ;
scatter x = q5_mean_theta_Tdu y = q5_mean_theta_Tms ;
run ;

