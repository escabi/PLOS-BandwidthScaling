%
%function [TRF,E]=trfalpha1model(beta,taxis)
%
%   FILE NAME       : TRF ALPHA 1 MODEL
%   DESCRIPTION     : Temporal receptive field model. This routine is
%                     similar to trfalphamodel.m but uses a 1-segment instead of a
%                     2-segment alpha function to definte the temporal recceptive
%                     field envelope
%
%   beta            : Parameter vector
%                     beta(1): Delay: Response latency (msec)
%                     beta(2): Fm: Character temporal modualtion frequency (Hz)
%                     beta(3): BW: Temporal modulation bandwidth (Hz)
%                     beta(4): Phase (recommended pi/4 for on-off neurons)
%                     beta(5): Peak time domain amplitude
%   taxis            : Time axis (msec) 
%
%RETURNED VARIABLES
%
%   TRF             : Model temporal receptive field (TRF)
%   E               : Envelope
%
% (C) F. He, M.A. Escabi, Feb/2018 (edit 03/21, 08/22)
%
function [TRF, E]=trfalpha1model(beta,taxis)

%Model Parameters
delay   = beta(1);
Fm      = beta(2);
bw      = beta(3);
tau     = 2*sqrt(sqrt(2)-1)/2/pi/bw*1000;% msec 
phase   = beta(4);
A       = beta(5);

K = 0;  %DC offset for alpha function model
E = alphafxn1model([delay,tau,1,K],taxis);
TRF = E.*cos(2*pi*Fm*(taxis-delay)/1000-phase); 
TRF = A.*TRF./max(TRF);

end
