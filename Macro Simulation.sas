/*############################################# Macro for Joint-Frailty-Simulation #############################################
################################################################################################################################*/

%macro SimulateJointFrailty(output,nrep,n,seed,frailtydist,theta,gamma,scalerec,shaperec,scaleterm,shapeterm,
							FU,censprob,beta1,beta2,path);

/* 	INPUT:
	output:		 	    Name of the output-dataset with the simulated values
						(The dataset contains nrep stacked single datasets, each with n patients)
						In addition, &output is the prefix of the summary-datasets describing the simulated data
	nrep:				Number of simulated single datasets
	n:					Number of subjects per dataset
	seed:				Seed for reproducible simulation 
	frailtydist:		Frailty distribution: gamma or lognormal
	theta:			    Variance of the frailty (mean is always 1) 
	gamma:				Exponent-parameter for the terminal event frailty 
	scalerec:			Scale-parameter of the Weibull-distribution for the recurrent events
	shaperec:			Shape-parameter of the Weibull-distribution for the recurrent events
	scaleterm:			Scale-parameter of the Weibull-distribution for the terminal event
	shapeterm:			Shape-parameter of the Weibull-distribution for the terminal event
	FU:					Maximum-length of the Follow-up (administrative censoring at that time)
	censprob:			Cumulative probability of random censoring (uniform distribution) in the interval [0,FU)
	beta1:				Effect of the covariate x (Bin[1,0.5]-distributed) on the recurrent event rate
	beta2:				Effect of the covariate x (Bin[1,0.5]-distributed) on the terminal event rate
	path:               path where to store the output-datasets (i.e. the dataset with simulated data and the 
						2 summary-datasets) in csv-format. If no output-datasets should be stored in csv-format, 
						please specify path = none  

	OUTPUT:
	- dataset &output containing the simulated stacked sub-datasets
	- dataset &output._summary1 containing the following statistics: Min,Mean,Med,Max of recurrent, terminal and censoring
	  event numbers over all simulated datasets
	- dataset &output._summary2 containing the following statistics: Min,Mean,Med,Max of subject-numbers with 0,1,2,3,4,>=5 
	  recurrent events over all simulated datsets
	If &path ^= none, the output-datasets are stored under &path in csv-format */

	%let Ntotal = %eval(&nrep*&n);

/* ------------------------------------------------------ Simulation --------------------------------------------------------------- */
	
 /* Start a data-step for data-generation */
	data &output;

	 /* Set seed */
		call streaminit(&seed);		

	 /* Start of generating patient-profiles */
		do i=1 to &Ntotal;

		 /* Introduce subjectid and sampleid */
			sampleid = ceil(i/&n);	 /* sampleid = dataset-number */
			subjectid = mod(i,&n);   /* subject-number within a dataset*/
			if subjectid = 0 then subjectid = &n;

		 /* Patient-specific Frailty */
			%if &theta = 0 %then							         
				%do;
				frailty = 1;
				%end;
			%else 													 
				%do;
				%if &frailtydist = lognormal %then 
					%do; 
					frailty = (1/sqrt(1 + &theta))*(rand("LOGNORMAL")**sqrt(log(1+&theta)));
					%end;
				%if &frailtydist= gamma %then 
					%do;
					frailty=&theta*rand("GAMMA",1/&theta);
					%end;
				%end;
		
		 /* Subject-specific covariate: Binary treatment indicator (Bin[1,0.5]-distributed) */	
			x = rand("BERNOULLI",0.5);      								

		 /* Subject-specific scale-parameters of the Weibull-distributions for
			recurrent and terminal events (accounting for frailty and covariates) */
			subjectscalerec  = &scalerec*frailty*exp(&beta1*x);  				
			subjectscaleterm = &scaleterm*(frailty**&gamma)*exp(&beta2*x);
		
		 /* Terminal event time */
			timeterm = rand("WEIBULL", &shapeterm, (1/subjectscaleterm)**(1/&shapeterm)); 

		/* Censoring time */
			%if &censprob = 0 %then
				%do;
				timecens = &FU;
				%end;
			%else
				%do;
				/* sample from a uniform distribution whose distribution function has value 1 at time (1/&censprob)*&FU,
				   i.e. the distribution function has slope &censprob/FU and hence its value at time &FU is &censprob.
				   Take minimum: This way the censoring distribution has &censprob mass on [0,&FU] and (1-&censprob) mass 
				   on &FU */
				timecens =  min((1/&censprob)*&FU*rand("UNIFORM"),&FU);
				%end;

		 /* Follow-up time */
			timeobserved = min(timeterm,timecens);						      			 							

		 /* Initialize Within-patient-loop */
			timestart = 0;			/* initialize starting time of the first at-risk-interval as 0 */
			timestop = 0;			/* initialize stopping time of the first at-risk-interval as 0 (will be overwritten later) */
			eventindicator = 1;		/* initialize eventindicator of the first at-risk-interval as 1 (will be overwritten later)
									   Remark: The eventindicator shows why the actual risk interval ended 
									   (0 for censoring, 1 for recurrent event, 2 for terminal event) */

		 /* Start of within-subject-loop: 
			Generate (probably multiple) event-rows for each subject in each dataset.
			Stop producing new recurrent event times if the last drawn recurrent 
			event time is larger than the follow-up-time, i.e. if the eventindicator
			of the last subject-row is 0 or 2 */	
			do until (eventindicator ne 1);																			

				timestart = timestop;	/* beginning of the new risk interval is the end of the old risk interval */																			
				A = rand("UNIFORM");
				u = -log(A);
				/* Time between last and new event time */
				intertime = ((u + subjectscalerec*(timestart**&shaperec))/subjectscalerec)**(1/&shaperec) - timestart;	
				/* New recurrent event time */
				timestop = timestart + intertime;	
				
				/* If new recurrent event is before end of follow-up, then a new recurrent event time is simulated 
				   (because eventindicator is initiated with 1) 
				   If new recurrent event is after end of follow-up, then  ... */	
				if (timestop > timeobserved) then 									
					do;
					timestop = timeobserved;		/* ... timestop is the end of follow-up */											
					if (timeterm <= timecens) then                                         
						eventindicator = 2;			/* ... eventindicator=2 if follow-up ends due to a terminal event */											
					else 
						eventindicator = 0;			/* ... eventindicator=0 if follow-up ends due to censoring */
					end;

				output;		/* Output dataset-row now */

			end; 	/* End of [do until (eventindicator ne 1)]-within-subject-loop */

		end;    /* Stop generating subject-profiles, i.e. the [do i=1 to &Ntotal]-loop */

	 /* Finally only sampleid, subjectid, frailty, x, timestart, timestop and eventindicator in dataset */
		keep sampleid subjectid  frailty x timestart timestop eventindicator;	

	run; /* End of the data-step for data-generation */


	/* ------------------------------------------------------ Make summary tables --------------------------------------------------------------- */

 /* Summary-Table 1: Distribution of recurrent, terminal, censoring event numbers within the nrep datasets */
	data &output._sort;
		set &output;
	run;

	proc sort data=&output._sort;
		by sampleid x;
	run;

	proc freq data=&output._sort noprint;
		by sampleid x;
		tables eventindicator / out=&output._summary1;
	run;

	proc means data=&output._summary1 noprint;
		var count;
		class eventindicator x;
		output out=&output._summary1(drop=_type_ _freq_) min=min mean=mean median=median max=max;
	run;

	data &output._summary1;
		set &output._summary1;
		label min = "min" mean = "mean" median = "median" max = "max";
		if (~missing(eventindicator)) and (~missing(x));
	run;


 /* Summary-Table 2: Distribution of patient numbers with 1,2,3,4,>=5 recurrent events within the nrep datasets */
	proc freq data=&output._sort noprint;
		by sampleid x;
		tables subjectid / out=&output._summary2;
	run;

	data &output._summary2;
		set &output._summary2(rename=(count=count1));
		count1=count1-1; /* because last line is always censoring or terminal event time */
		if count1 >=5 then Recevents = '>=5';
		else Recevents=count1;
		drop percent;
	run;

	proc freq data=&output._summary2 noprint;
		by sampleid x;
		tables Recevents / out=&output._summary2;
	run;

	proc means data=&output._summary2 noprint;
		var count;
		class Recevents x;
		output out=&output._summary2(drop=_type_ _freq_) min=min mean=mean median=median max=max;
	run;

	data &output._summary2;
		set &output._summary2;
		label min = "min" mean = "mean" median = "median" max = "max";
		if (~missing(Recevents)) and (~missing(x));
	run;

	proc datasets;
	   delete &output._sort;
	quit;

/* ------------------------------------------------------ Export --------------------------------------------------------------- */

	%if &path ^= none %then
		%do;
		proc export data=&output 
   			outfile="&path.\&output..csv"
   			dbms=dlm
   			replace;
			delimiter=';';
		run;
		proc export data=&output._summary1 
   			outfile="&path.\&output._summary1.csv"
   			dbms=dlm
   			replace;
			delimiter=';';
		run;
		proc export data=&output._summary2 
   			outfile="&path.\&output._summary2.csv"
   			dbms=dlm
   			replace;
			delimiter=';';
		run;
		%end;

%mend SimulateJointFrailty;
