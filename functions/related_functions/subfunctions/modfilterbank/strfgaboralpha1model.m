%
%function [STRFm]=strfgaboralpha1model(beta,input)
%
%   FILE NAME       : STRF GABOR ALPHA1 MODEL
%   DESCRIPTION     : Separable STRF model. The spectral receptive field
%                     is modeled as a gabor function (gaussian multiplied with cosine)
%                     while the temporal receptive field is modeled as the 
%                     product of a 1-segment alpha fucntion and a cosine.
%
%                     Note that the STRFs from this model differ from strfgaboralpha1.m and
%                     strfgaboralpha2.m since they are generated using 1-segmenet alpha
%                     functions as opposed to 2-segment alpha functions.alpha. In those instances,
%                     1 and 2 designate the number of separable STRFs used in the model.
%
%   beta            : STRF parameter vector
%                     beta(1): Response latency (msec)
%                     beta(2): Fm: Character Modulation Frequency (Hz)
%                     beta(3): Temporal modulation bandwidth (Hz)
%                     beta(4): Temporal phase (0-2*pi)/defualt pi/4
%                     beta(5): Best octave frequency, xo
%                     beta(6): Best spectral modulation frequency (cycles/octaves)
%                     beta(7): Spectral Modulation Bandwidth (cycles/octaves)
%                     beta(8): Spectral phase (0-2*pi)
%                     beta(9): Amplitude
%   input.taxis     : Time axis (msec)
%   input.X         : Octave frequency axis (octaves)
%
%RETURNED VARIABLES
%
%   STRFm           : Separable STRF model
%
% (C) F. He, M.A. Escabi, Feb/2018 (Last Edit 05/21)
%
function [STRFm]=strfgaboralpha1model(beta,input)

%Parameter vectors
betat=[beta(1:4) 1];    % Parameters for trfalpha1model
betas=[beta(5:8) 1];    % Parameters for srfgabormodel
betas(2)=1/betas(2);    % Convert spectral modulation bandwidth (cycles/oct) to Receptive Field Bandwidth (oct)

%Generating the separable STRF model
[TRF]=trfalpha1model(betat,input.taxis);
[SRF]=srfgabormodel(betas,input.X);
STRFm = beta(9) * SRF'*TRF;
end
