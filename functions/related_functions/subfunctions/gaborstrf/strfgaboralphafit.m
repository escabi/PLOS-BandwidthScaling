%
%function [GAModel]=strfgaboralphafit(STRFs,STRF,taxis,faxis,N,theta,betaLB,betaUB,L)
%
%   FILE NAME       : STRF GABOR ALPHA FIT
%   DESCRIPTION     : STRF model opimization routine. The STRF is fitted
%                     using a least squares minimization procedure. The
%                     temporal and spectral receptive fields are fitted
%                     with a gabor alpaha function.
%
%   STRFs           : Significant STRF
%   STRF            : Original STRF
%   taxis           : Time axis (sec)
%   faxis           : Frequency Axis (Hz)
%   N               : Maximum number of iterations (global search)
%   theta           : Target error percentage
%   betaLB          : Lower bound for model parameters (Optional) - assume
%                     1st order model
%   betaUB          : Upper bound for model parameters (Optional) - assumes
%                     1st order model
%   L               : Maximum number of function evaluations for LSQCURVEFIT
%                     (Default = 500)
%
%RETURNED VARIABLES
%
%   GAModel           : Data structure containing the following model
%                      results
%   .taxis           : Time axis (msec)
%   .faxis           : Frequency Axis
%   .STRFm1          : First order STRF model
%   .STRFm2          : Second order STRF model
%   .beta1           : STRF parameter vector
%        PARAMETERS FOR FIRST STRF COMPONENTS
%                     beta(1): Response latency (msec)
%                     beta(2): Rise time constant (msec)
%                     beta(3): Decay time constant (msec)
%                     beta(4): Best temporal modulation frequency (Hz)
%                     beta(5): Temporal phase (0-2*pi)
%                     beta(6): Best octave frequency, xo
%                     beta(7): Best spectral modulation frequency (cycles/octaves)
%                     beta(8): Gaussian spectral bandwidth (octaves)
%                     beta(9): Spectral phase (0-2*pi)
%                     beta(10): Peak Amplitude
%   .beta2           : The second order nonseparable model also includes 
%       PARAMETERS FOR SECOND STRF COMPONENT ALSO INCLUDE
%                     beta(11): Response latency (msec)
%                     beta(12): Rise time constant (msec)
%                     beta(13): Decay time constant (msec)
%                     beta(14): Best temporal modulation frequency (Hz)
%                     beta(15): Temporal phase (0-2*pi)
%                     beta(16): Best octave frequency, xo
%                     beta(17): Best spectral modulation frequency (cycles/octaves)
%                     beta(18): Gaussian spectral bandwidth (octaves)
%                     beta(19): Spectral phase (0-2*pi)
%                     beta(20): Peak Amplitude
%   .Cov1           : Parameter covariance matrix for first order model
%   .Cov2           : Parameter covariance matrix for second order model
%   .FI1            : Fisher information matrix for first order model
%   .FI2            : Fisher information matrix for second order model
%   .P1             : Ratio test for first order model: The model is 
%                     significant if P(10)>1.96
%   .P2             : Ratio test for second order model: The model is
%                     significant if P(10)>1.96 & P(20)>1.96
%   .MSE1           : Normalized mean square error for model 1 - noise
%                     variance removed (e.g., Zheng & Escabi 2013)
%   .MSE2           : Normalized mean square error for model 2 - noise
%                     variance removed (e.g., Zheng & Escabi 2013)
%   .SI1            : Similarity index betwen STRFm1 & STRFs       
%   .SI2            : Similarity index between STRFm2 & STRFs
%   .N1             : Number of optimization itterations for Model 1
%   .N2             : Number of optimization itterations for Model 2
%   .Order          : Model order determined by significance test
%
% (C) Created June 2023 by Jenna P. Blain(Modified from Monty A. Escabi, October 2006. 
%     Revised by Chen Chen, Mar. 2007 (Edit MAE
%     Feb. 2017; MAE May 2021). Based on original code by Anqui Qiu 2003.
%
function [GAModel]=strfgaboralphafit(STRFs,STRF,taxis,faxis,N,theta,betaLB,betaUB,L)

%Paramater Lower and Upper bounds
if nargin<7 | isempty(betaLB)
    betaLB=[0  0  0   0    0 0 0 0 0    0                    ];
end
if nargin<8 | isempty(betaUB)
    betaUB=[50 20 50 500 2*pi 8 4 8 2*pi 5*max(max(abs(STRF)))];
end
if nargin<9 | isempty(L)
    L=500;
end

%Estimating initial STRF parameters
[RFParam]=strfparam(taxis-min(taxis),faxis,STRFs,500,4);


beta(1)=RFParam.PeakDelayP/2;
beta(2)=RFParam.PeakDelayP/2;
beta(3)=RFParam.Duration10-RFParam.PeakDelay/2;
beta(4)=RFParam.bTMF;
beta(5)=pi/4;
beta(6)=RFParam.BF;
beta(7)=mean(abs(RFParam.BestRD));
if beta(7)<1/betaUB(8)
   beta(8)=1;
else
   beta(8)=1/beta(7);
end
beta(9)=0;
beta(10)=max(max(STRF));

%Maximum number of evaluations
lsqOpts = optimoptions('lsqcurvefit','MaxFunctionEvaluations', L);

%Fitting separable STRF model (to estimate initial parameters with
%significant STRF)
input.taxis=1000*(taxis-min(taxis));
input.X=log2(faxis/faxis(1));
[beta10,RESNORM1,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,J1]=lsqcurvefit('strfgaboralpha1',beta,input,STRFs,betaLB,betaUB,lsqOpts);

%Fitting separable STRF model to real data - beta from prior step is used
%for initial parameters
[beta1,RESNORM1,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,J1]=lsqcurvefit('strfgaboralpha1',beta10,input,STRF,betaLB,betaUB,lsqOpts);
STRFm1=strfgaboralpha1(beta1,input);
index1=find(STRFs~=0);
index2=find(STRFs==0);
MSE1=( var(STRFm1(index1)-STRFs(index1)) - var(STRF(index2)) ) / (var(STRFs(index1))) * 100;

%Fitting nonseparable STRF model to real data and finding normalized mean squared erorr
betaLB2=[betaLB betaLB];
betaUB2=[betaUB betaUB];
beta20=[beta1 beta1];
[beta2,RESNORM2,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,J2]=lsqcurvefit('strfgaboralpha2',beta20,input,STRF,betaLB2,betaUB2,lsqOpts);
STRFm2=strfgaboralpha2(beta2,input);
index1=find(STRFs~=0);
index2=find(STRFs==0);
MSE2=( var(STRFm2(index1)-STRFs(index1)) - var(STRF(index2)) ) / (var(STRFs(index1))) * 100;

%If Separable model error is above threshold, randomize some key parameters and reoptimize
n1=1;
while MSE1>theta && n1<N

       beta=beta1.*(1+(rand(size(beta10))-0.5)*0.25);       %Add 25% error randomly
       beta(4)=rand*betaUB(4);                              %Best temporal modulation frequency - fully randomize within search space
       beta(5)=2*pi*rand;                                   %Temporal phase - fully randomize within search space
       beta(7)=rand*betaUB(7);                              %Ripple density - fully randomize within search space
       if beta(7)<1/betaUB(8)
          beta(8)=1;
       else
          beta(8)=1/beta(7).*(1+(rand-0.5)*0.25);           %Spectral bandwidth - change so that it is consistent with initial ripple density.
       end
       beta(9)=2*pi*rand;                                   %Spectral phase - fully ranodmize within search space
      
       [beta,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,J]=lsqcurvefit('strfgaboralpha1',beta,input,STRF,betaLB,betaUB,lsqOpts);
      
           if RESNORM<RESNORM1
              J1=J;
              beta1=beta;
              STRFm1=strfgaboralpha1(beta1,input);
              MSE1=( var(STRFm1(index1)-STRFs(index1)) - var(STRF(index2)) ) / (var(STRFs(index1))) * 100;
              RESNORM1=RESNORM;
           end   
        
       n1=n1+1;
           
end

%If Nonseparable model error is above threshold, randomize some key
%parameters and reoptimize
n2=1;
while MSE2>theta && n2<N
    

     beta=beta2.*(1+(rand(size(beta20))-0.5)*0.25);         %Add 25% error randomly
       beta(4)=rand*betaUB2(4);                              %Best temporal modulation frequency - fully randomize within search space
       beta(5)=2*pi*rand;                                   %Temporal phase - fully randomize within search space
       beta(7)=rand*betaUB2(7);                              %Ripple density - fully randomize within search space
       if beta(7)<1/betaUB2(8)
          beta(8)=1;
       else
          beta(8)=1/beta(7).*(1+(rand-0.5)*0.25);           %Spectral bandwidth - change so that it is consistent with initial ripple density.
       end
       beta(9)=2*pi*rand;                                   %Spectral phase - fully ranodmize within search space
  
       beta(14)=rand*betaUB2(14);                             %Best temporal modulation frequency - fully randomize within search space
       beta(15)=2*pi*rand;                                   %Temporal phase - fully randomize within search space
       beta(17)=rand*betaUB2(17);                             %Ripple density - fully randomize within search space
       if beta(17)<1/betaUB2(18)
          beta(18)=1;
       else
          beta(18)=1/beta(17).*(1+(rand-0.5)*0.25);           %Spectral bandwidth - change so that it is consistent with initial ripple density.
       end
       beta(19)=2*pi*rand;                                   %Spectral phase - fully ranodmize within search space
 
      
       [beta,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,J]=lsqcurvefit('strfgaboralpha2',beta,input,STRF,betaLB2,betaUB2,lsqOpts);
      
           if RESNORM<RESNORM2
              J2=J;
              beta2=beta;
              STRFm2=strfgaboralpha2(beta2,input);
              MSE2=( var(STRFm2(index1)-STRFs(index1)) - var(STRF(index2)) ) / (var(STRFs(index1))) * 100;
              RESNORM2=RESNORM;
           end   
        
       n2=n2+1;
           
end 

%Covariance Matrix
Cov1=full(inv(J1'*J1));         %Covariance Matrix
Cov2=full(inv(J2'*J2));         %Covariance Matrix
index=find(STRFs==0);
Var=var(STRF(index));           %Noise variance estimate
FI1=Cov1*Var;                   %Fisher Information Matrix
FI2=Cov2*Var;                   %Fisher Information Matrix
P1=beta1./sqrt(diag(FI1)');     %Ratio Test
P2=beta2./sqrt(diag(FI2)');     %Ratio Test

% Calculating similarity index
[RSTRF1] = strfcorr(STRFs,STRFm1,taxis,faxis);
[RSTRF2] = strfcorr(STRFs,STRFm2,taxis,faxis);

%Creating Data Structure
GAModel.taxis=taxis;
GAModel.faxis=faxis;
GAModel.STRFs=STRFs;
GAModel.STRF=STRF;
GAModel.STRFm1=STRFm1;
GAModel.STRFm2=STRFm2;
GAModel.beta1=beta1;
GAModel.beta2=beta2;
GAModel.Cov1=Cov1;
GAModel.Cov2=Cov2;
GAModel.FI1=FI1;
GAModel.FI2=FI2;
GAModel.P1=P1;
GAModel.P2=P2;
GAModel.MSE1=MSE1;
GAModel.MSE2=MSE2;
GAModel.SI1=RSTRF1.SI;
GAModel.SI2=RSTRF2.SI;
GAModel.N1=n1;
GAModel.N2=n2;

%Significance test to determine order
if P2(20)>1.96 & P2(10)>1.96
    GAModel.Order=2;
elseif P1(10)>1.96;
    GAModel.Order=1;
else
    GAModel.Order=0;
end