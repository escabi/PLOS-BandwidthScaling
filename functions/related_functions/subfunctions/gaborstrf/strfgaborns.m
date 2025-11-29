%
%function [STRFm]=strfgaborsm(beta,input)
%
%   FILE NAME       : STRF GABOR NS
%   DESCRIPTION     : Non Separable STRF model. The receoptive field
%                     envelope is modeled as a two dimensional Gaussian and 
%                     the carrier component uses a nonseparable sinusoid.
%
%   beta            : STRF parameter vector
%                     beta(1): Peak delay (msec)
%                     beta(2): Gaussian temporal duration (msec)
%                     beta(3): Best temporal modulation frequency (Hz)
%                     beta(4): Best octave frequency, xo
%                     beta(5): Gaussian spectral bandwidth (octaves)
%                     beta(6): Best spectral modulation frequency (octaves)
%                     beta(7): Spectro-Temporal Phase (0-2*pi)
%                     beta(8): Peak Amplitude
%   input.taxis     : Time axis (msec)
%   input.X         : Octave frequency axis (octaves)
%
%RETURNED VARIABLES
%
%   STRFm           : Speraable STRF model
%
% (C) Monty A. Escabi, Aug 2017
%
function [STRFm]=strfgaborsm(beta,input)

%Time and Frequency Axis
taxis=input.taxis;
X=input.X;
X=(repmat(X',1,length(taxis)));
taxis=repmat(taxis,size(X,1),1)/1000;

%Parameters
Del=beta(1)/1000;
D=beta(2)/1000;
Fm=beta(3);
X0=beta(4);
BW=beta(5);
RD=beta(6);
P=beta(7);
A=beta(8);


%Spectro-Tempioral Gabor Model
STRFm=A*exp(-(2*(X-X0)/BW).^2).*exp(-(2*(taxis-Del)/D).^2).*cos(2*pi*RD*(X-X0)+2*pi*Fm*(taxis-Del)+P);
