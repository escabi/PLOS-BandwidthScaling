%
%function [MidData]=midbrainogram(CochData,beta,Sflag)
%
%   FILE NAME   : MIDBRAIN OGRAM
%   DESCRIPTION : Auditory midbrain model sound representation obtained 
%                 by applying an octave spaced modualtion filterbank on the
%                 a cochleogram representation. Takes the peripheral
%                 auditory model decomposition (cochleogram) as input. The
%                 input cochleogram is then filtered by a spectroteporal
%                 modulation filterbank that mimics the auditory midbrain
%                 decomposition. The resulting multi-dimensional output
%                 represents the sound in time, frequency, temporal
%                 moduatlion frequency and spectral modualtion frequency.
%
%                 If the primary goal is to estiamte the modulations in
%                 the sound through a midbrain model (as in He et al), it
%                 is suggested that the nonlinearity used in the cochleogram 
%                 represenation is set to NLType='hil'. If on the other hand
%                 the goal is to implement a more realistic model of the 
%                 auditory periphery and midbrain, it is suggested that the 
%                 nonlinearity used in the cochleogram represenation is set 
%                 to NLType='rect' which rectifies the input and models the 
%                 fine struture phase locking for low frequencies in the 
%                 cochlea.
%
%                 For high frequency cochlear channels (e.g., >1kHz), the
%                 'rect' and 'hil' nonlinearities produce similar outputs. 
%                 In both instances, they extract the lowpass filtered
%                 envelopes. Note that the carriers are discarded for the high
%                 frequency channels because the synaptic lowpass filter of
%                 the hair cell membrane (modeled by the synaptic lowpass
%                 filter Fm cutoff) limits synchornization to Fm Hz (which
%                 is typically set between 0.75 and 1.0 kHz). By compariosn
%                 for low frequencies (e.g., <1kHZ), the two nonlinearties 
%                 behave differently and the 'rect' is a arguably the more 
%                 realistic option. The 'hil' option extracts only the evelope 
%                 of the low frequency channels and it ignores the fact 
%                 that low frequency channels also syncrhonize to the carrier
%                 information ("fine structure"), which the hilbert 
%                 transform discards. On the other hand, the 'rect' option 
%                 allows for both the carrier information (i.e., the "fine 
%                 strcuture") and the envelopes of these low frequency channels 
%                 to be encoded, which is more realistic from a periphery 
%                 and midbrain coding perspective.
%
%   CochData    : Output cochleogram - see 'cochleogram.m' for details
%   beta        : Modulation filterbank parameters vector (see ModFiltBankParamGaborAlpha1ns for details)
%                     beta(1): Qt quality factor for temporal modualtion filters   (recommended: 1)
%                     beta(2): Qs quality factor for spectral modulation filters  (recommended: 1)
%                     beta(3): Fml - lower temporal modulation frequency (Hz)
%                     beta(4): Fmu - upper temporal modulation frequency (Hz)
%                     beta(5): RDl - lower spectral modulation frequency (cycles/oct)
%                     beta(6): RDu - upper spectral modulation frequency (cycles/oct)
%                     beta(7): Dt  - temporal modulation filter spacing (octave)
%                     beta(8): Ds  - spectral modulation filter spacing (octave)
%   Sflag       : Determines which spectrogram represenation to decompose:
%                 'dB', 'dBR' or 'Lin' (Optional, Default == 'Lin'). 'dBR' 
%                 is the dB cochleogram with a compressibe nonlinearity that rescales
%                 the amplitudes from -1 to 1. Check out COCHLEOGRAM for details.
%                 
%RETURNED VARIABLES
%
%   MidData : Data structure containing midbrainogram results
%             .Sm                   : Midbrain model output represenation.
%                                     Multidimensional output NxLxMxN where
%
%                                     N=number of frequency channels
%                                     L=number of time samples
%                                     M=number of spectral modulation channels
%                                     N=number of temporal modulation channels
%                           
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
% (C) Monty A. Escabi, Feb 2022 (Edit 10/22)
%
function [MidData]=midbrainogram(CochData,beta,Sflag)

%Input Arguments
if nargin<3
    Sflag='Lin';
end

%Selecting cochleogram represenation - either dB, rescaled dB, or Linear amplitude
%Selecting additional parameters
if strcmp(Sflag,'dB')
    S=CochData.SdB;         %dB amplitude cochleogram
elseif strcmp(Sflag,'dBR')
    S=CochData.SdBR;        %Rescaled dB amplitude cochleogram
else
    S=CochData.S;           %Linear amplitude cochleogram
end
Fst=1/CochData.taxis(2);
X=log2(CochData.faxis/CochData.faxis(1));
Fss=1/X(2);

%Generating STRF Modualtion Filter Bank
[ModBankParam] = ModFiltBankParamGaborAlpha1ns(beta);               %Generating filter parameters
[STRFBank] = STRFModFiltBankGaborAlpha1ns(ModBankParam,[Fst Fss]);  %Generating STRF filters

%Initializing Convolution Matrix
Lt=size(S,2);
Ls=size(S,1);
MaxNt=-9999;
for k=1:size(STRFBank.F,1)
    for l=1:size(STRFBank.F,2)
            MaxNt=max(size(STRFBank.F(k,l).H,2),MaxNt); %Finding Maximum STRF duration
    end
end
Sm=zeros(Ls,Lt+MaxNt-1,size(STRFBank.F,1),size(STRFBank.F,2));

%Decomposiong Cochleogram into Spectrotemporal Modualtion Bands
clc
for k=1:size(STRFBank.F,1)
    for l=1:size(STRFBank.F,2)

            %Convolving, truncating spectral edges, and storing 
            Y=conv2(S,STRFBank.F(k,l).H);     
            Ns=(size(STRFBank.F(k,l).H,1)-1)/2;     %1/2 STRF spectral response size
            Nt=size(STRFBank.F(k,l).H,2);           %STRF temporal response size
            Sm(:,1:Lt+Nt-1,k,l)=Y(Ns+1:end-Ns,:);   %Truncating spectral edges
         
    end
        %Output Display
        clc
        disp(['Computing STRF Modulation Decomposition: ' num2str(k/size(STRFBank.F,1)*100,3) ' %'])
end

%Computing Midbrain Modulation Power Spectrum
mMPS=squeeze(mean(Sm.^2,2));                        %MAE, 10/22

%Storing Results in Data Structure
MidData.Sm=Sm;
MidData.mMPS=mMPS;
MidData.taxis=(0:size(Sm,2)-1)/Fst;
MidData.faxis=CochData.faxis;
MidData.Fmaxis=[ModBankParam.F(1,:).Fm];
MidData.RDaxis=[ModBankParam.F(:,1).RD];
MidData.Param.ModBankParam=ModBankParam;
MidData.Param.CochParam=CochData.Param;
MidData.Param.beta=beta;
MidData.Param.Sflag=Sflag;