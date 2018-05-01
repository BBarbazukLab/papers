libname trago '/ufrc/barbazuk/lboat/Old_World_New_BED/sas_data';

/* for parents */
%macro make_sbys_parents (reads, alt) ;

/* Transpose dataset so reps in column  */
    proc sort data=trago.ase_4_bayes_&reads._reads_tdu_&alt. ;
        by commonID rep;
        run;

    proc transpose data=trago.ase_4_bayes_&reads._reads_tdu_&alt. out=flip_total prefix=total_count_rep;
        by commonID ;
        var total_count;
        id rep;
        run;
    proc transpose data=trago.ase_4_bayes_&reads._reads_tdu_&alt. out=flip_tdu prefix=tdu_count_rep;
        by commonID ;
        var tdu_count;
        id rep;
        run;
    proc transpose data=trago.ase_4_bayes_&reads._reads_tdu_&alt. out=flip_&alt. prefix=&alt._count_rep;
        by commonID ;
        var &alt._count;
        id rep;
        run;

data trago.ase_bayes_&reads._rds_tdu_&alt._sbys ;
retain commonID ;
merge flip_total flip_tdu flip_&alt. ;
by commonID ;
drop _NAME_ ;
run ;
%mend ;

%make_sbys_parents (tdu, tpo);
%make_sbys_parents (tpo, tpo);
%make_sbys_parents (tdu, tpr);
%make_sbys_parents (tpr, tpr);
