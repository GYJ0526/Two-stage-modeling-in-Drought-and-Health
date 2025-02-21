libname NOAA 'C:\Users\yeongjin.gwon\OneDrive - University of Nebraska Medical Center\01 Academics\01 Research\04 Spatial Modeling\02 NOAA\1. Manuscript\5. Upper Midwest\3. Data';

/*******************************************/
/* Overall analysis for general population */
/*******************************************/
proc import out=NOAA.test
  datafile= "C:\Users\yeongjin.gwon\OneDrive - University of Nebraska Medical Center\01 Academics\01 Research\04 Spatial Modeling\02 NOAA\1. Manuscript\5. Upper Midwest\3. Data\UpperMidwest.csv"
  dbms=csv replace;
  getnames=yes;
  datarow=2;
run;

/* Descriptive statistics */
proc means data=NOAA.test min p25 p50 p75 max mean maxdec=3;
  class State Year;
  var TD CD RD;
run;

proc means data=NOAA.test min p25 p50 p75 max mean sum maxdec=3;
  class county;
  var race;
run;

proc means data=NOAA.test min p25 p50 p75 max mean maxdec=3;
  class State;
  var EDDIB06 EDDIB12;
run;

proc means data=NOAA.test min p25 p50 p75 max mean maxdec=3;
  class Year;
  var Dmod Dsev SPEIM06 SPEIM12 SPEIS06 SPEIS12 EDDIM06 EDDIM12 EDDIS06 EDDIS12;
run;

/*******************************/
/* Bayesian Two-stage Modeling */
/*******************************/
%macro BayesTwostage(OUTCOME,D1,D2,Drought);
data test;
  set NOAA.test;
  idc+1;
  by State County Year Month;
  if first.County then idc=1;
  Lpop=log(Pop);
*  if Pop>=12500 then output;
run;

/* Stage1. Quasi-Poisson regression model */
ods exclude all;
proc glimmix data=test;
  by county;
  class Year;
*  effect spl=spline(idc/ basis=bspline details);
  effect spl=spline(idc/ degree=6 knotmethod=percentile(10));
  model &OUTCOME=&D1 &D2 spl Year TempA TempA*TempA / link=log offset=Lpop solution;
  _variance_ = _mu_;
  random _residual_;
  ods output ParameterEstimates=_para_est_USDM ConvergenceStatus=ConvergeStatus;
run;
ods exclude none;

data OUT;
  merge _para_est_USDM ConvergeStatus;
  by County;
  if Status=0 then output;
run;

proc sort data=OUT out=OUT_USDM_M;
  by Effect;
run;

data Moderate;
  set OUT_USDM_M(where=(Effect=&Drought));
  SE2=0;
  if Estimate<=-1.5 or Estimate>=1.0 then delete;
  if StdErr>=1.0 then delete;
  Est=round(Estimate,0.001);
  SE=round(StdErr,0.001);
  SE2=SE**2;
  if SE=.  then delete;
  keep County Est SE SE2 Probt;
run;

/* Stage2. Bayesian linear regression */
proc mcmc data=Moderate outpost=postoutB nmc=50000 nbi=10000 nthin=5 seed=1234 propcov=quanew DIC maxtune=30;
 
  parms alpha phi;
  prior alpha~normal(0,var=100);
  prior phi~normal(0,var=5);
  random muK~normal(0,var=1) subject=County;
  xalpha=alpha+muK*exp(phi);
  model Est~normal(xalpha,var=SE2);

run;

data NOAA.twostage;
  set postoutB;
  Group="TS";
  RR=exp(alpha);
run;

proc univariate data=NOAA.twostage;
  var RR;
  output OUT=PostSummary MEAN=mean STD=sd PCTLPTS=2.5 97.5 PCTLPRE=P;
run;

proc print data=PostSummary;
run;

/* Compute Posterior probability for the Risk ratio */
data Postprob;
  set NOAA.twostage;
  keep RR P05 P10 P15 P20 P30;
  P05=0;P10=0;P15=0;P20=0;P30=0;
  if RR>=1.05 then P05=1;
  if RR>=1.10 then P10=1;
  if RR>=1.15 then P15=1;
  if RR>=1.20 then P20=1;
  if RR>=1.30 then P30=1;
run;

proc report data=Postprob nowd center split='|' headline headskip 
  style(report)=[font_face='Times New Roman' font_size=11pt] 
  style(column)=[font_face='Times New Roman' font_size=11pt]
  style(lines) =[font_face='Times New Roman' font_size=11pt just=c font_weight=bold foreground=black background=GRAYB8];

  columns P05 P10 P15 P20 P30;
  define  P05 / "RR>=5% "           mean format=6.4;
  define  P10 / "RR>=10% "           mean format=6.4;
  define  P15 / "RR>=15% "           mean format=6.4;
  define  P20 / "RR>=20% "           mean format=6.4;
  define  P30 / "RR>=30%"           mean format=6.4;
run;

%mend;

%BayesTwostage(OUTCOME=RD,D1=EDDIM06,D2=EDDIS06,Drought="EDDIS06");
%BayesTwostage(OUTCOME=RD,D1=EDDIM12,D2=EDDIS12,Drought="EDDIS12");
%BayesTwostage(OUTCOME=RD,D1=EDDIM01,D2=EDDIS01,Drought="EDDIS01");








/*********************************************/
/* Overall graphical assessment scatter plot */
/*********************************************/
data NOAA.full;
  set NOAA.histout NOAA.currout NOAA.Power NOAA.twostage;
  Variance=exp(phi);
run;

proc sgplot data=NOAA.full;
    scatter x=RR y=phi / group=Group;
run;

proc sgplot data=NOAA.full;
   vbox RR / group=Group;
   keylegend / title="Method";
run;

proc export data=NOAA.full
  outfile="c:/test/PosteriorSample.csv"
  dbms=csv replace;
run;




proc import out=NOAA.testrace
  datafile= "C:\Users\yeongjin.gwon\OneDrive - University of Nebraska Medical Center\01 Academics\01 Research\04 Spatial Modeling\02 NOAA\1. Manuscript\5. Upper Midwest\3. Data\UpperMidwest_Race.csv"
  dbms=csv replace;
  getnames=yes;
  datarow=2;
run;

/* Descriptive statistics */
proc sort data=NOAA.testrace out=testrace;
  by State County race;
run;

proc means data=testrace max mean sum maxdec=1;
  class race;
  var Pop;
run;

proc means data=testrace max mean sum maxdec=1;
  class State race;
  var Pop;
run;
