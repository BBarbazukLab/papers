/*Test for power-law properties*/
proc import datafile='/ufrc/barbazuk/jobrant/SPINY/SAS_check/aco_mus_network_data.csv' 
            out=aco_mus_network_data
            dbms=csv 
            replace; 
        getnames=yes; 
run; 

ods graphics on;

/* Get the frequency values */
proc freq data=work.aco_mus_network_data noprint;
    table Degree / out=degree_freq;
run;

/*Print Procedure -- with head 10*/
proc print data=degree_freq(obs=10); 
run;

proc sgplot data=work.aco_mus_network_data;
    histogram Degree / legendlabel="Node Degree Histogram";
run;

data raw_degree_freq;
    set degree_freq;
    Key = Degree;
    Value = COUNT;
run;

proc export data=raw_degree_freq
            outfile="/ufrc/barbazuk/jobrant/SPINY/SAS_check/node_freq_table.csv"
            dbms=csv 
            replace; 
run; 

/* Get log frequency and degree */
data log_degree_freq;
    set degree_freq;
    logDegree = log10(Degree);
    logFreq = log10(COUNT);
run;

/*Table of Contents Procedure*/
proc contents data=log_degree_freq short; 
run;

/*Print Procedure -- with head 10*/
proc print data=log_degree_freq(obs=10); 
run; 

/*Display Univariate Summary Statistics -- logDegree*/
proc univariate data=log_degree_freq;
        var logDegree; 
run; 

/* Perform linear regression for power-law */
proc reg data=log_degree_freq plots=diagnostics(unpack); 
    model logFreq=logDegree; 
run;

/* Plot power-law */
proc sgplot data=log_degree_freq;
    reg X=logDegree Y=logFreq;
run;



