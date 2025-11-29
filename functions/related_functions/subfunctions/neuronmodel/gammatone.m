%
%function [X]=gammatone(beta,time)
%
%   FILE NAME       : GAMMA TONE 
%   DESCRIPTION     : Geneartes a gammatone
%
%   beta            : Parameter vector containing
%           beta(1) : a, amplitude
%           beta(2) : N, gammatone order
%           beta(3) : b, bandwidth (Hz)
%           beta(4) : fc, center frequency (Hz)
%           beta(5) : d, delay (sec)
%
%   time            : Time axis vector (sec)
%
%RETURNED VARIABLES
%
%	X               : Gammatone function
%
% (C) Escabi & Clonan, Aug 2024
%
function [X]=gammatone(beta,time)

%REFERENCES:    Van Eemeerseel & Peeters, Acoustic Research Letters 2003
%               de Boer 1975 (First developed the GTF)

%Input Parameters
a=beta(1);
N=beta(2);
b=beta(3);
fc=beta(4);
d=beta(5);

%Generating Impulse Response
X=zeros(size(time));
i=find(time>=d);
X(i)=a*(time(i)-d).^(N-1).*exp(-2*pi*b*(time(i)-d)).*sin(2*pi*fc*(time(i)-d));