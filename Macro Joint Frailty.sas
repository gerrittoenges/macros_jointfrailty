
/*############################################# Macro for Joint-Frailty-Analysis ###############################################
################################################################################################################################*/

%macro JointFrailty(input,output,sampleidvar,subjectidvar,timevar,eventindicatorvar,linpredrec,linpredterm,frailtydist,
					methodgamma,hazards,startval,optimstartval,quad,quadpoints,path);

/*	INPUT:
	input: 				name of the input-dataset containing the simulated event-times 
	output: 			prefix for all output dataset-names 
	sampleidvar: 		In case of the stacked datasets of a simulation: Name of the variable in the input-dataset 
						that contains the sample-ID. If only one single dataset is analysed and no sample-ID-variable 
						exists: please specify sampleidvar = none
	subjectidvar: 		Name of the variable in the input-dataset that contains the subject-ID
	timevar: 			Name of the variable in the input-dataset that contains the event-times
	eventindicatorvar: 	Name of the variable in the input-dataset that indicates the type of event 
						(0=censoring, 1=recurrent, 2=terminal) 	
	linpredrec: 		Linear predictor for the recurrent events; e.g. in case of the two covariates 
						x1 and x2 you have to specify linpredrec = beta11*x1 + beta12*x2. Categorical 
						covariates must be split up into dummy-variables in advance
	linpredterm: 		Linear predictor for the terminal event; e.g. in case of the two covariates x1 and x2 
						you have to	specify linpredterm = beta21*x1 + beta22*x2. Categorical covariates must be 
						split up into dummy-variables in advance
	frailtydist: 		Frailty distribution (gamma or lognormal)
	methodgamma: 		If frailtydist = lognormal then set methodgamma = none. Otherwise, if frailtydist = gamma, 
						you have to choose the method for dealing with the non-random effect (i.e. lr or pit)
	hazards: 			Either constant, weibull or piecewise (10 pieces)
	startval: 			name of the (input-)dataset that contains the starting values for the likelihood-maximization. 
	optimstartval: 		Should starting values be optimized by first fitting a simplified model in NLMIXED without
						frailty term (true or false)?
	quad: 				Which numerical quadrature-technique should be used?
						(ad=adaptive gaussian quadrature, noad=non-adaptive gaussian quadrature)
	quadpoints: 		Either "auto" for letting NLMIXED determine the number of quadrature points itself or
			    		a number which specifies how many quadrature points should be used 
	path: 				path where to store the output-datasets in csv-format. If no output-datasets should be stored 
						in csv-format, please specify path = none

	OUTPUT:
	- dataset &output._est containing the parameter estimates of all sub-datasets
	- dataset &output._conv containing the convergence-status (0 = convergence or 1 = failed convergence) of all 
	  sub-datasets
	- dataset &output._time containing the processing time for the nlmixed-procedure
	If &path ^= none, the output-datasets are stored under &path in csv-format */



/* ---------------------------  Dataset-Preparation in case of analysis with piecewise constant hazards ----------------------- */

	%if &hazards = piecewise %then
		%do;

		/* Create a sample-ID-variable (even in case of &sampleidvar = none) */
		%if &sampleidvar = none %then 
			%do;
			data &input._for_JointFrailty;
				sampleID_temp_ = 1;
				set &input;
			run;
			%end;
		%else 
			%do;
			data &input._for_JointFrailty;
				set &input;
				sampleID_temp_ = &sampleidvar;
			run;
			%end;

		/* Determine the largest Follow-up time */
		proc univariate data=&input._for_JointFrailty noprint;
			by sampleID_temp_;
			var &timevar;
			output out=maxval max=max;
		run;

		/* Obtain quantiles for recurrent events */
		proc univariate data=&input._for_JointFrailty(where=(&eventindicatorvar=1)) noprint;
			by sampleID_temp_;
			var &timevar;
			output out=quant_r pctlpts= 10 20 30 40 50 60 70 80 90 pctlpre=qr;
		run;
		data quant_r;
			retain sampleID_temp_ qr0 qr10 qr20 qr30 qr40 qr50 qr60 qr70 qr80 qr90 max;
			merge quant_r maxval;
			by sampleID_temp_;
			rename  qr10=qr1 qr20=qr2 qr30=qr3 qr40=qr4 qr50=qr5 qr60=qr6 qr70=qr7 qr80=qr8 qr90=qr9 max=qr10;
			qr0=0;
			label qr10=' ' qr20=' ' qr30=' ' qr40=' ' qr50=' ' qr60=' ' qr70=' ' qr80=' ' qr90=' ' max=' ';
		run;	/* qr0 is 0, qr10 is the largest follow-up time, qr1-qr9 are the 10%,...,90%-quantiles of recurrent events */

		/* Obtain quantiles for death */
		proc univariate data=&input._for_JointFrailty(where=(&eventindicatorvar=2)) noprint;
			by sampleID_temp_;
			var &timevar;
			output out=quant_d pctlpts= 10 20 30 40 50 60 70 80 90 pctlpre=qd;
		run;
		data quant_d;
			retain sampleID_temp_ qd0 qd10 qd20 qd30 qd40 qd50 qd60 qd70 qd80 qd90 max;
			merge quant_d maxval;
			by sampleID_temp_;
			rename  qd10=qd1 qd20=qd2 qd30=qd3 qd40=qd4 qd50=qd5 qd60=qd6 qd70=qd7 qd80=qd8 qd90=qd9 max=qd10;
			qd0=0;
			label qd10=' ' qd20=' ' qd30=' ' qd40=' ' qd50=' ' qd60=' ' qd70=' ' qd80=' ' qd90=' ' max=' ';
		run;	/* qd0 is 0, qd10 is the largest follow-up time, qd1-qd9 are the 10%,...,90%-quantiles of death events */
			
		/* Calculate the duration in each interval and the indicator of event in each interval */
		data &input._for_JointFrailty;
			merge &input._for_JointFrailty quant_r quant_d;
			by sampleID_temp_;
		run;

		data &input._for_JointFrailty;
			set &input._for_JointFrailty; 
			drop sampleID_temp_;
			array quant_r {*} qr0-qr10;
			array quant_d {*} qd0-qd10;
			array dur_r {*} dur_r1-dur_r10;
			array dur_d {*} dur_d1-dur_d10;			    
			array event_r {*} event_r1-event_r10;
			array event_d {*} event_d1-event_d10;
			    
			do i=1 to 10; 		/* Initialization with 0's */
				dur_r{i}=0;
				dur_d{i}=0;
				event_r{i}=0;
				event_d{i}=0;
			end;
			    
			if &eventindicatorvar=1 then 
				do; /* For line with recurrent event: likelihood-contribution is just the recurrent hazard */
			   		do i=2 to 11;
			  			if &timevar <= (1/1000000)*ceil(1000000*quant_r{i}) then 
			            	do;
			                	event_r{i-1}=1; /* only one event_r-Variable is 1; dur_r-, event_d-, dur_d-Variables are all 0 */
			                	i=11;
							end;
			       	end;
				end;
			else 
				do; /* For line with terminal event or censoring
					   likelihood-contribution censoring: cumulative recurrent hazard, cumulative terminal hazard 
					   likelihood-contribution terminal event: cumulative recurrent hazard, cumulative terminal hazard, terminal hazard */
			  		do i=2 to 11;
						if &timevar <= (1/1000000)*ceil(1000000*quant_d{i}) then 
							do;
			                	event_d{i-1}=(&eventindicatorvar=2); /* censoring: all event_d-Variables are 0; terminal event: only one event_d-Variable is 1 */
			                	dur_d{i-1}=&timevar-quant_d{i-1};    /* censoring and terminal event: Only the dur_d-Variable of the interval that contains the (censoring/terminal) event time gets a duration-entry */
			                	i=11;								
			            	end;
			            else dur_d{i-1}=quant_d{i}-quant_d{i-1}; /* censoring and terminal event: All dur_d-Variables of intervals preceding the interval that contains the (censoring/terminal) event time get a duration-entry */
					end;
					do i=2 to 11;
						if &timevar <= (1/1000000)*ceil(1000000*quant_r{i}) then              
							do;
			                	dur_r{i-1}=&timevar-quant_r{i-1}; /* censoring and terminal event: Only the dur_r-Variable of the interval that contains the (censoring/terminal) event time gets a duration-entry */
			                	i=11;
			            	end;
						else dur_r{i-1}=quant_r{i}-quant_r{i-1}; /* censoring and terminal event: All dur_r-Variables of intervals preceding the interval that contains the (censoring/terminal) event time get a duration-entry */
					end;
				end;
		run;

		proc datasets;
	   		delete quant_r quant_d maxval;
		quit;

		%end;

	/* In case of constant or weibull baseline-hazards we can directly use the input-dataset */
	%else
		%do;
		data &input._for_JointFrailty;
			set &input;
		run;
		%end;

/* -------------------------------------------- Optimization of starting values ------------------------------------------------- */

	%if &optimstartval = true %then
		%do;

		data startval_temp1;
		 	set &startval;
			where parameter ^= "theta" and parameter ^= "gamma";
		run;

		data startval_temp2;
		 	set &startval;
			where parameter = "theta" or parameter = "gamma";
		run;

		proc nlmixed data=&input._for_JointFrailty;

			/* Starting values */
			parms / data=startval_temp1;
	
			/* Separate analyses for each sample */
			%if &sampleidvar ^= none %then
				%do;
				by &sampleidvar;
				%end;

			/* Parameter bounds */
			%if &hazards = piecewise %then 
				%do;
				bounds r01 r02 r03 r04 r05 r06 r07 r08 r09 r10 h01 h02 h03 h04 h05 h06
				h07 h08 h09 h10 >= 0;
				%end; 
			%if &hazards = constant %then
				%do;
				bounds r01 h01 >= 0;
				%end;
			%if &hazards = weibull %then
				%do;
				bounds scalerec shaperec scaleterm shapeterm >= 0;
				%end;

			/* baseline hazards */
			%if &hazards = piecewise %then 
				%do;
				base_haz_r = r01*event_r1 + r02*event_r2 + r03*event_r3 + r04*event_r4 +
				r05*event_r5 + r06*event_r6 + r07*event_r7 + r08*event_r8 
				+ r09*event_r9 + r10*event_r10;
				cum_base_haz_r = r01*dur_r1 + r02*dur_r2 + r03*dur_r3 + r04*dur_r4 + 
				r05*dur_r5 + r06*dur_r6 + r07*dur_r7 + r08*dur_r8 + r09*dur_r9 + 
				r10*dur_r10;
				base_haz_d = h01*event_d1 + h02*event_d2 + h03*event_d3 + h04*event_d4 + 
				h05*event_d5 + h06*event_d6 + h07*event_d7 + h08*event_d8 + 
				h09*event_d9 + h10*event_d10;
				cum_base_haz_d = h01*dur_d1 + h02*dur_d2 + h03*dur_d3 +
				h04*dur_d4 + h05*dur_d5 + h06*dur_d6 + h07*dur_d7 +
				h08*dur_d8 + h09*dur_d9 + h10*dur_d10;
				%end;
			%if &hazards = constant %then 
				%do;
				base_haz_r = r01;											
				cum_base_haz_r = r01 * &timevar;  								
				base_haz_d = h01;												
				cum_base_haz_d = h01 * &timevar;     						   
				%end;
			%if &hazards = weibull %then 
				%do;
				base_haz_r = scalerec*shaperec*&timevar**(shaperec-1);		
				cum_base_haz_r = scalerec*&timevar**shaperec; 				
				base_haz_d = scaleterm*shapeterm*&timevar**(shapeterm-1);	
				cum_base_haz_d = scaleterm*&timevar**shapeterm;   			    						    
				%end;

			/* Linear predictors (recurrent, terminal) WITHOUT random effect */
			linpred1= &linpredrec;     
			linpred2= &linpredterm;            

			loglik1=-exp(linpred1) * cum_base_haz_r;
			loglik2=-exp(linpred2) * cum_base_haz_d;
				
			/* Log-likelihood-contribution of each dataset-row */
			if &eventindicatorvar = 1 then loglik=log(base_haz_r) + linpred1;  
			if &eventindicatorvar = 2 then loglik=loglik1 + log(base_haz_d) + linpred2 + loglik2;
			if &eventindicatorvar = 0 then loglik=loglik1 + loglik2;

			/* Model statement */
			model &timevar ~ general(loglik);

			/* Select output-tables */
			ods output ParameterEstimates=optim_startval1;	 
 
		run; /* End of the nlmixed-procedure */

		%if &sampleidvar ^= none %then 
			%do;
			data optim_startval1;
				length parameter $20;
				set optim_startval1;
				keep &sampleidvar parameter estimate;
				call symputx('sample_number',&sampleidvar);
			run;

			data _null_; 
				call symputx('list',repeat('startval_temp2 ',&sample_number - 1)); 
			run;
			data optim_startval2; 
				set &list;
				&sampleidvar + 1; 
			run;
			data optim_startval2; 
				length parameter $20;
				set optim_startval2;
				&sampleidvar = ceil(&sampleidvar/2);
			run;

			data optim_startval;
				length parameter $20;
				set optim_startval1 optim_startval2;
				by &sampleidvar;
			run;
			%end;
		%else 
			%do;
			data optim_startval1;
				length parameter $20;
				set optim_startval1;
				keep parameter estimate;
			run;

			data optim_startval2;
				length parameter $20;
				set startval_temp2;
			run;

			data optim_startval;
				length parameter $20;
				set optim_startval1 optim_startval2;
			run;

			%end;


		run;
		
		%end;



/* ----------------------------------------- Joint Frailty model fit with NLMIXED-Procedure -------------------------------------------- */

	%let timer_start = %sysfunc(datetime()); /* Start timer */

	/* Which type of numerical quadrature, how many quadrature-points */
	%if &quad=ad %then
		%do;
		%if &quadpoints^=auto %then
			%do;
			proc nlmixed data=&input._for_JointFrailty qpoints=&quadpoints itdetails;
			%end;
		%else
			%do;
			proc nlmixed data=&input._for_JointFrailty itdetails;
			%end;		 
		%end;
	 %if &quad=noad %then
		%do;
		%if &quadpoints^=auto %then
			%do;
			proc nlmixed data=&input._for_JointFrailty qpoints=&quadpoints itdetails noad;
			%end;
		%else
			%do;
			proc nlmixed data=&input._for_JointFrailty itdetails noad;
			%end; 
		%end;
				
		/* Starting values */
	 %if &optimstartval = true %then
	 	%do;
		parms / bydata data=optim_startval;
		%end;
	%else
		%do;
		parms / data=&startval;
		%end;

		/* Separate analyses for each sample */
		%if &sampleidvar ^= none %then
			%do;
			by &sampleidvar;
			%end;
				
		/* Parameter bounds */
		%if &hazards = piecewise %then 
			%do;
			bounds r01 r02 r03 r04 r05 r06 r07 r08 r09 r10 h01 h02 h03 h04 h05 h06
			h07 h08 h09 h10 theta >= 0;
			%end; 
		%if &hazards = constant %then
			%do;
			bounds r01 h01 theta >= 0;
			%end;
		%if &hazards = weibull %then
			%do;
			bounds scalerec shaperec scaleterm shapeterm theta >= 0;
			%end;

			    
		/* baseline hazards */
		%if &hazards = piecewise %then 
			%do;
			base_haz_r = r01*event_r1 + r02*event_r2 + r03*event_r3 + r04*event_r4 +
			r05*event_r5 + r06*event_r6 + r07*event_r7 + r08*event_r8 
			+ r09*event_r9 + r10*event_r10;
			cum_base_haz_r = r01*dur_r1 + r02*dur_r2 + r03*dur_r3 + r04*dur_r4 + 
			r05*dur_r5 + r06*dur_r6 + r07*dur_r7 + r08*dur_r8 + r09*dur_r9 + 
			r10*dur_r10;
			base_haz_d = h01*event_d1 + h02*event_d2 + h03*event_d3 + h04*event_d4 + 
			h05*event_d5 + h06*event_d6 + h07*event_d7 + h08*event_d8 + 
			h09*event_d9 + h10*event_d10;
			cum_base_haz_d = h01*dur_d1 + h02*dur_d2 + h03*dur_d3 +
			h04*dur_d4 + h05*dur_d5 + h06*dur_d6 + h07*dur_d7 +
			h08*dur_d8 + h09*dur_d9 + h10*dur_d10;
			%end;
		%if &hazards = constant %then 
			%do;
			base_haz_r = r01;											
			cum_base_haz_r = r01 * &timevar;  								
			base_haz_d = h01;												
			cum_base_haz_d = h01 * &timevar;     						   
			%end;
		%if &hazards = weibull %then 
			%do;
			base_haz_r = scalerec*shaperec*&timevar**(shaperec-1);		
			cum_base_haz_r = scalerec*&timevar**shaperec; 				
			base_haz_d = scaleterm*shapeterm*&timevar**(shapeterm-1);	
			cum_base_haz_d = scaleterm*&timevar**shapeterm;   			    						    
			%end;


		/* If Gamma-Frailty*/
		%if &frailtydist = gamma %then 
			%do;		

			/* If Likelihood-reformulation-method */
			%if &methodgamma = lr %then
				%do;
				/* log-likelihood-appending summands */
				/* 1. logofloggamma:
					  log(f(nu)) where f(nu) is the density of log(X) and X is Gamma-distributed with 
					  E(X)=1 and Var(X)=theta 
					  Attention: Do not use 0 as starting value for theta (dividion through 0!) */
				logofloggamma = (1/theta)*nu - (1/theta)*exp(nu) - lgamma(1/theta) - (1/theta)*log(theta); 
				/* 2. logofstandardnormal:
				      log(f(nu)) where f(nu) is the standard normal density */ 
				logofstandardnormal = -(nu**2)/2;                                                                

				/* Linear predictors (recurrent, terminal) */
				linpred1= &linpredrec + nu;             
			  	linpred2= &linpredterm + gamma * nu;      

				loglik1=-exp(linpred1) * cum_base_haz_r;
				loglik2=-exp(linpred2) * cum_base_haz_d;

				/* Log-likelihood-contribution of each dataset-row */
				if &eventindicatorvar = 1 then loglik=log(base_haz_r) + linpred1;  
				if &eventindicatorvar = 2 then loglik=loglik1 + log(base_haz_d) + linpred2 + loglik2 + logofloggamma - logofstandardnormal; 
				if &eventindicatorvar = 0 then loglik=loglik1 + loglik2 + logofloggamma - logofstandardnormal;

				/* Model statement */
				model &timevar ~ general(loglik);

				/* Random statement */
				random nu ~ normal(0,1) subject=&subjectidvar;
				%end;


			/* If methodgamma = pit */
			%if &methodgamma = pit %then
				%do;
				p = cdf('NORMAL',nu,0,1);
				if p>0.999999 then p=0.999999;
				g = quantile('GAMMA',p,1/theta,theta);
				/* g is Gamma-distributed with mean 1 and variance theta */

				/* Linear predictors (recurrent, terminal) */
				linpred1= &linpredrec + log(g);             
			  	linpred2= &linpredterm + gamma * log(g);     

				loglik1=-exp(linpred1) * cum_base_haz_r;
				loglik2=-exp(linpred2) * cum_base_haz_d;

				/* Log-likelihood-contribution of each dataset-row */
				if &eventindicatorvar = 1 then loglik=log(base_haz_r) + linpred1;  
				if &eventindicatorvar = 2 then loglik=loglik1 + log(base_haz_d)+linpred2 + loglik2; 
				if &eventindicatorvar = 0 then loglik=loglik1 + loglik2;

				/* Model statement */
				model &timevar ~ general(loglik);

				/* Random statement */
				random nu ~ normal(0,1) subject=&subjectidvar;
				%end;

			%end;
			

		/* If Lognormal-Frailty */
		%if &frailtydist = lognormal %then
			%do;															   
			/* Linear predictors (recurrent, terminal) */
			linpred1= &linpredrec + nu;             
		  	linpred2= &linpredterm + gamma * nu;      

			loglik1=-exp(linpred1) * cum_base_haz_r;
			loglik2=-exp(linpred2) * cum_base_haz_d;

			/* Log-likelihood-contribution of each dataset-row */
			if &eventindicatorvar = 1 then loglik=log(base_haz_r) + linpred1;  
			if &eventindicatorvar = 2 then loglik=loglik1 + log(base_haz_d) + linpred2 + loglik2; 
			if &eventindicatorvar = 0 then loglik=loglik1 + loglik2;

			/* Model statement */
			model &timevar ~ general(loglik);

			/* Random statement */
			random nu ~ normal(log(1/sqrt(1+theta)),log(1+theta)) subject=&subjectidvar;
			%end;


		/* Select output-tables */
		ods output ParameterEstimates=&output._est;
		ods output ConvergenceStatus=&output._conv;	 
 
	run; /* End of the nlmixed-procedure */
	

/*-----------------------------------------------  Export output-datasets  -----------------------------------------------------*/

	/* Delete the _for_JointFrailty-dataset */
	proc datasets;
   		delete &input._for_JointFrailty;
	quit;

	/* Delete the datasets in case of optimstartval = true */
	%if &optimstartval = true %then
		%do;
		proc datasets;
   			delete optim_startval optim_startval1 optim_startval2 startval_temp1 startval_temp2;
		quit;
		%end;

	/* Stop timer */ 
	data &output._time;
  		duration_in_min = (datetime() - &timer_start)/60;
	run;
	
	/* Export output-datasets in csv-format if desired */
	%if &path ^= none %then
		%do;
		/* Export the nlmixed-processing-time */
		proc export data=&output._time 
	   		outfile="&path.\&output._time.csv"
	   		dbms=dlm
	   		replace;
			delimiter=';';
		run;
		/* Export the the whole simulation results */
		proc export data=&output._est
	   		outfile="&path.\&output._est.csv"
	   		dbms=dlm
	   		replace;
			delimiter=';';
		run;
		proc export data=&output._conv
	   		outfile="&path.\&output._conv.csv"
	   		dbms=dlm
	   		replace;
			delimiter=';';
		run;
		%end;


%mend JointFrailty;
