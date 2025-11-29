%
%function [SRF,E]=srfgabormodel(beta,x)
%
%   FILE NAME       : SRF GABOR MODEL
%   DESCRIPTION     : Spectral receptive field gabor model
%
%   beta            : Gabor parameter vector
%                     beta(1): Best octave frequency, xo
%                     beta(2): Best spectral modulation frequency (octaves)
%                     beta(3): Gaussian spectral bandwidth (octaves)
%                     beta(4): Spectral phase (0-2*pi) 
%                     beta(5): Peak amplitude 
%   x               : Octave frequency axis (octaves)
%
%RETURNED VARIABLES
%
%   SRF             : Model spectral receptive field (SRF)
%   E               : TRF model envelope
%
% (C) Monty A. Escabi, October 2006 (Edit 5/21 MAE; Edit 8/21 MAE)
%
function [SRF,E]=srfgabormodel(beta,x);

E=exp(-(2*(x-beta(1))/beta(3)).^2);
SRF=beta(5)*E.*cos(2*pi*beta(2)*(x-beta(1))+beta(4));