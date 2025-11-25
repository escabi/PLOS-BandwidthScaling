%
%function [MidRateData]=midbrainrateogram(MidData,fc,OF)
%
%	FILE NAME 	: MIDBRAIN RATE O GRAM
%	DESCRIPTION : Midbrain rate represenation of a sound. Takes auditory
%	              midbrain model data as input (see midbrainogram.m). Each
%	              output of the midbrain model is rectified and lowpass
%	              filtered with a exponential filter with cutoff fc that 
%                 mimimics a cell membrane integration for a rate coding
%                 strategy.
%
%   MidData     : Output midbrainogram.m
%   fc          : Lowpass cutoff frequency (Hz)
%   OF          : Oversampling factor (Default == 1)
%                 
%RETURNED VARIABLES
%
%   MidRateData - Data structure containing ctxogram results
%             .Sm                   : Midbrain rate model output represenation.
%                                     Multidimensional output NxLxMxN where
%
%                                     N=number of frequency channels
%                                     L=number of time samples
%                                     M=number of spectral modulation channels
%                                     N=number of temporal modulation channels
%                           
%             .midMPS               : Mid modulation power spectrum.
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
function [MidRateData] = midbrainrateogram(MidData,fc,OF)

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
midMPS=squeeze(mean(Y.^2,2));                   %MAE, 10/22

%Appending Data to Structure
MidRateData=MidData;                                %Compying midbrain data structure since it contains relevant information
MidRateData.Sm=Y;                                   %Replacing output with Mid Model Output
MidRateData.midMPS=midMPS;                          %Mid Modulation Power Spectrum
MidRateData.taxis=(0:size(Y,2)-1)/(Fs/DF);          %Resampled time axis
MidRateData.Param.MidRateParam.fc=fc;               %MidRate fc parameter
MidRateData.Param.MidRateParam.OF=OF;               %MidRate Tau parameters

end