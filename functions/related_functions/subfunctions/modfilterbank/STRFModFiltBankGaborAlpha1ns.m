%
%function [STRFFiltBank] = STRFModFiltBankGaborAlpha1ns(ModBankParam,Fsd)
%
%   FILE NAME       : STRF Modualtion Filerbank Gabor Alpha1 Nonseparable
%
%   DESCRIPTION     : Generates a filterbank of non-separable Gabor
%                     (spectral) and 1st order alpha function (temporal)
%                     model STRF receptive fields
%
%   ModBankParam.F(j,i) - Data structure Matrix containg Modulation 
%                         Filterbank STRF parameters. The Data structure is
%                         generated using ModFiltBankParamGaborAlpha1ns.m.
%                         It contains the following fields:
%
%   .Delay              : Response latency (msec)
%   .Fm                 : Character Modulation Frequency (Hz) (3dB)
%   .BWt                : Temporal modulation bandwidth (Hz)
%   .Bof                : Best octave frequency, xo (default = 0)
%   .RD                 : Best spectral modulation frequency (cycles/octaves)
%   .BWs                : Spectral Modulation Bandwidth (cycles/octaves) (sqrt(2)/(2*sigma))
%   .P                  : Spectro-temporal phase (0-2*pi) (default pi/4)
%   .Amp                : Amplitude (excludes carrier)
%   .Beta               : Everything above in sequential order
%   
%   ModBankParam.Param - Modulation Filter Bank Parameter Data Structure
%
%   .Qt                 : Qt quality factor for temporal modualtion filters
%   .Qs                 : Qs quality factor for spectral modualtion filters
%   .Fml                : Lower temporal modulation frequency (Hz)
%   .Fmu                : Upper temporal modulation frequency (Hz)
%   .RDl                : Lower spectral modulation frequency (cycles/oct
%   .RDu                : Upper spectral modulation frequency (cycles/oct)
%   .Dt                 : Temporal modulation filter spacing (octave)
%   .Ds                 : Spectral modulation filter spacing (octave)
%   .beta               : Original input vector containing filterbank parameters
%
%RETURNED VARIABLES
%
%   STRFFiltBank - Data structure containg modualtion filterbank STRFs. It
%                  contains the following fields:
%
%   .F(j,i).H           : Spectro-Temporal Receptive Field Model Filters
%   .F(j,i).input.taxis : Time Axis Vector (ms)
%   .F(j,i).input.X     : Frequency Axis Vector (octave)
%   .F(j,i).beta        : STRF parameters
%
% (C) M.A. Escabi, Oct. 2021
%
function [STRFFiltBank] = STRFModFiltBankGaborAlpha1ns(ModBankParam,Fsd)

%Tempporal and Spectral Threshold - used to select STRF matrix size
Tt=7;           %Number of time constants 
Ts=4;           %Number of standard deviations

%Temporal and Spectral Sampling Rates
Fs_t=Fsd(1);
Fs_s=Fsd(2);

%Generate STRF Filterbank
[Ns,Nt] = size(ModBankParam.F);
for i = 1:Nt
    for j = 1:Ns
        %Compute STRF Time constant and Bandwidth
        Beta = ModBankParam.F(j,i).Beta;                    %STRF parameters
        tau = 2*sqrt(sqrt(2)-1)/2/pi/Beta(3);               %Convert Temporal Modualtion Bandwdith to STRF Time constant
        bw = 4*sqrt(log(2))/pi/Beta(6);                     %Convert Spectral Modulation Bandwidth to STRF bandwidth
        
        %Computing Nonseparable STRF model
        input.taxis=(0:ceil(tau*Tt*Fs_t))/Fs_t*1000;
        input.X=(-ceil(bw*Ts*Fs_s):ceil(bw*Ts*Fs_s))/Fs_s;
        H = strfgaboralpha1modelns(Beta,input);
         
        %Normalize the STRF to achieve unit passband gain
        As=2/bw/sqrt(pi);                                   %Spectral Normalization Gain
        At=pi/exp(1)/sqrt(sqrt(2)-1)*Beta(3);               %Temporal Normalization Gain
        H=H*As*At/Fs_t/Fs_s;                                %Normalize by Temporal and Spectral sampling rates in order to get correct gain
        
        %Assing STRF to STRF Filter Bank structure
        STRFFiltBank.F(j,i).H=H;                    
        STRFFiltBank.F(j,i).input=input;
        STRFFiltBank.F(j,i).beta=Beta;
    end
end