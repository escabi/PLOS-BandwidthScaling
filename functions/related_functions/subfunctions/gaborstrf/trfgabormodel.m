%
%function [TRF,E]=trfgabormodel(beta,taxis)
%
%   FILE NAME       : TRF GABOR MODEL
%   DESCRIPTION     : Temporal receptive field gabor model
%
%   beta            : Gabor parameter vector
%                     beta(1): Peak delay (msec)
%                     beta(2): Gaussian temporal duration (msec)
%                     beta(3): Best temporal modulation frequency (Hz)
%                     beta(4): Temporal phase (0-2*pi)
%                     beta(5): Time warping coefficient
%                     beta(6): Peak amplitude
%   taxis           : Time axis (msec)
%
%RETURNED VARIABLES
%
%   SRF             : Model temporal receptive field (TRF)
%   E               : SRF model envelope
%
% (C) Monty A. Escabi, October 2006 (edit Jul 2018)
%
function [TRF,E]=trfgabormodel(beta,taxis);

t=2*atan(beta(5)*taxis/1000)-beta(1)/1000;
E=exp(-(2*t/(beta(2)/1000)).^2);                %Return envelope (MAE, 2018)
TRF=beta(6)*E.*cos(2*pi*beta(3)*t+beta(4));