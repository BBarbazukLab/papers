
libname trago '/ufrc/barbazuk/lboat/Old_World_New_BED/sas_data';

/* Create flag_analyze -- cast has 4 sets */

%macro flag_analyze (reads, alt) ;

/* cast reads */
    * A commonID is called analyzable if it has ASE information for at least one rep.;
    proc sort data=trago.ase_4_bayes_&reads._reads_tdu_&alt. ;  
        by commonID;
        run;

    /* Removed noprint from end of means for checks */
    proc means data=trago.ase_4_bayes_&reads._reads_tdu_&alt. noprint;
        by commonID;
        output out=sums sum(ASE_count)=sums;
        run;

    * needed to add _freq_ requirement to ensure that there are 3 reps.;
    data flag_analyze;
        set sums; 
        if sums > 0 and _freq_ ge 3 then flag_analyze = 1; else flag_analyze = 0;
        keep commonID flag_analyze;
        run;

    proc print data=flag_analyze (obs=10);
    run;

/* Make permanent dataset */
    data trago.bayes_flag_&reads._reads_tdu_&alt. ;
        set flag_analyze;
        run;

%mend ;

%flag_analyze (Tm, tpo) ;
%flag_analyze (Tms, tpr) ;
