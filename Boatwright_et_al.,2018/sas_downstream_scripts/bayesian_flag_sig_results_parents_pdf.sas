
libname trago "/ufrc/barbazuk/lboat/Old_World_New_BED/sas_data";

/* parents */
%macro flag (reads, ref1, ref2) ; ;

data trago.bayes_flag_sig_&reads._for_&ref1._&ref2;
    set trago.Emp_&reads._for_&ref1._&ref2 ;
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

%flag (tdu, tdu, tpo) ;
%flag (tpo, tdu, tpo) ;
%flag (tdu, tdu, tpr) ;
%flag (tpr, tdu, tpr) ;

/* plot parents */
%macro plots (reads, ref1, ref2) ;

ods pdf file="/ufrc/barbazuk/lboat/Old_World_New_BED/reports/bayes_flag_sig_&reads._for_&ref1._&ref2..pdf";

title "Distribution of mean (*ESC*){unicode theta} for Emp_&reads._for_&ref1._&ref2";
proc sgplot data =  trago.bayes_flag_sig_&reads._for_&ref1._&ref2. ;
histogram q4_mean_theta / fillattrs=graphdata1 legendlabel="q = 0.4";
density q4_mean_theta / lineattrs=graphdata1 legendlabel="q4 density";
histogram q5_mean_theta / fillattrs=graphdata2 legendlabel="q = 0.5";
density q5_mean_theta / lineattrs=graphdata2 legendlabel="q5 density";
histogram q6_mean_theta / fillattrs=graphdata3 legendlabel="q = 0.6";
density q6_mean_theta / lineattrs=graphdata3 legendlabel="q6 density";
xaxis label="mean (*ESC*){unicode theta}";
run ;

ods pdf close ;
%mend ;

%plots (tdu, tdu, tpo) ;
%plots (tpo, tdu, tpo) ;
%plots (tdu, tdu, tpr) ;
%plots (tpr, tdu, tpr) ;

