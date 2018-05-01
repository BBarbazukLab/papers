
libname trago "";


/* parents */
%macro flag (reads, ref1, ref2) ; ;

data trago.bayes_flag_sig_&reads._for_&ref1._&ref2;
    set Emp_&reads._for_&ref1._&ref2 ;
	if (q4_q025 ne " ") and (q4_q025 > 0.5 or q4_q975 < 0.5) then flag_q4_sig = 1; 
		else flag_q4_sig = 0;
	if (q5_q025 ne " ") and (q5_q025 > 0.5 or q5_q975 < 0.5) then flag_q5_sig = 1; 
		else flag_q5_sig = 0;
    if (q6_q025 ne " ") and (q6_q025 > 0.5 or q6_q975 < 0.5) then flag_q6_sig = 1; 
		else flag_q6_sig = 0;
	if flag_q4_sig = 1 and flag_q5_sig=1 and flag_q6_sig = 1 then flag_sig_&reads._&ref1._&ref2 =1 ;
		else flag_sig_&reads._&ref1._&ref2 =0;
	*keep commonID flag_: ;
	run ;


    run;
%mend ;

%flag (Tdu, Tdu, tpo) ;
%flag (Tdu, Tdu, tpr) ;

%flag (Tpo, Tdu, tpo) ;
%flag (Tpr, Tdu, tpr) ;



/* plot parents */
%macro plots (reads, ref1, ref2) ;
title "Distribution of mean_theta for Emp_&reads._for_&ref1._&ref2";
proc sgplot data = trago.bayes_flag_sig_&reads._for_&ref1._&ref2;
histogram q4_mean_theta / fillattrs=graphdata1 ;
density q4_mean_theta / lineattrs=graphdata1;
histogram q5_mean_theta / fillattrs=graphdata2 ;
density q5_mean_theta / lineattrs=graphdata2;
histogram q6_mean_theta / fillattrs=graphdata3 ;
density q6_mean_theta / lineattrs=graphdata3;
xaxis label="mean_theta";
run ;
%mend ;

%plots (Tdu, Tdu, tpo) ;
%plots (Tdu, Tdu, tpr) ;

%plots (Tpo, Tdu, tpo) ;
%plots (Tpr, Tdu, tpr) ;

