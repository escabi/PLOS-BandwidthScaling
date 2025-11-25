%
%function [CtxData]=ctxogram(MidData,fc,OF)
%
%	FILE NAME 	: CTX O GRAM
%	DESCRIPTION : Cortexogram represenation of a sound. Takes auditory
%	              midbrain model data as input (see midbrainogram.m). Each
%	              output of the midbrain model is rectified and lowpass
%	              filtered with a exponential filter with cutoff fc that 
%                 mimimics a cell membrane integration. 
%
%   MidData     : Output midbrainogram.m
%   fc          : Lowpass cutoff frequency (Hz)
%   OF          : Oversampling factor (Default == 1)
%                 
%RETURNED VARIABLES
%
%   CtxData - Data structure containing ctxogram results
%             .Sm                   : Ctx model output represenation.
%                                     Multidimensional output NxLxMxN where
%
%                                     N=number of frequency channels
%                                     L=number of time samples
%                                     M=number of spectral modulation channels
%                                     N=number of temporal modulation channels
%                           
%             .ctxMPS               : Ctx modulation power spectrum.
%             .taxis                : Time axis
%             .faxis                : Frequency axis
%             .RDAxis               : Ripple density axis (cycles/oct)
%             .FmAxis               : Temporal modulation axis (Hz)
%             .STRFBank             : Modualtion filterbank STRFs
%             .Param.ModBankParam   : Modualtion filterbank parameters
%             .Param.CochParam      : Cochleogram parameters
%             .Param.beta           : Input hyper-parameters used to generate
%                                     modulation filterbank
%             .Param.Sflag          : Flag used to deterime type of
%                                     cochleogram
%
% (C) Fengrong He / Monty A. Escabi, June 2022
%
function [CtxData] = ctxogram(MidData,fc,OF)

%Input Arguments
if nargin < 3
    OF=1;
end

%Extracting data
Sm = MidData.Sm;
Fs = 1/(MidData.taxis(2)-MidData.taxis(1));

%Generating Lowpass Exponential Filter - mimics cell membrane filter
H = expfilter(Fs,fc);

%Downsampling Factor 
DF=floor(Fs/(fc*2)/OF);

%Rectifying and Lowpass Filtering
for k = 1:size(Sm,1)
    for l = 1:size(Sm,3)
        for m = 1:size(Sm,4)
   
            %Rectifying
            X = squeeze(Sm(k,:,l,m));
            X = max(0,X);                       %Rectification 
        
            %Lowpass filtering and downsampling to desired sampling rate
            X = conv(X,H);                      %Lowpass Filtered Output
            Y(k,:,l,m) = X(1:DF:end);           %Downsampling by DF

        end
    end
end

%Computing Midbrain Modulation Power Spectrum
ctxMPS=squeeze(mean(Y.^2,2));                   %MAE, 10/22

%Appending Data to Structure
CtxData=MidData;                                %Compying midbrain data structure since it contains relevant information
CtxData.Sm=Y;                                   %Replacing output with Ctx Model Output
CtxData.ctxMPS=ctxMPS;                          %Ctx Modulation Power Spectrum
CtxData.taxis=(0:size(Y,2)-1)/(Fs/DF);          %Resampled time axis
CtxData.Param.CtxParam.fc=fc;                   %Ctx fc parameter
CtxData.Param.CtxParam.OF=OF;                   %Ctx Tau parameters

end