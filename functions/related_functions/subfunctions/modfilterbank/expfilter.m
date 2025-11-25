%
%function [H] = expfilter(Fs,fc,NT)
%
%       FILE NAME       : EXP FILTER
%       DESCRIPTION     : First order exponential filter impulse response
%                         
%	Fs		: Sampling rate
%	fc		: Lowpass cutoff frequency (Hz)
%	NT      : Number of time constants for filter duration (Default NT=5)
%
%OUTPUT SIGNAL
%
%	H       : Impulse response
%
% (C) Fengrong He / M. Escabi, June 2022
%
function [H] = expfilter(Fs,fc,NT)

%Input Args
if nargin<3
    NT = 5; 
end

%Generating Impulse Response
Tau = 1/2/pi/fc;
N = ceil(NT*Tau*Fs);
time = (0:N-1)/Fs;
H=1/Tau*exp(-time/Tau)/Fs;

end
