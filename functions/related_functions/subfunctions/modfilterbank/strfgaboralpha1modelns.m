%
%function [STRFm]=strfgaboralpha1modelns(beta,input)
%
%   FILE NAME       : STRF GABOR ALPHA1 MODEL NS
%   DESCRIPTION     : Non-separable STRF model. The spectral receptive field
%                     is modeled as a gabor function (gaussian multiplied with cosine)
%                     while the temporal receptive field is modeled as the
%                     product of a 1-segment alpha fucntion and a cosine. The envelope of
%                     this model is seaparable, however, the carrier component
%                     (phase) is nonseparable. This allows the receptive fields to be selective
%                     for + and - spectro-temporal modulations.
%
%                     Note that the STRFs from this model differ from strfgaboralpha1.m and
%                     strfgaboralpha2.m since they are generated using 1-segmenet alpha
%                     functions as opposed to 2-segment alpha functions.alpha. In those instances,
%                     1 and 2 designate the number of separable STRFs used in the model.
%
%                     This model is similar to strfgaboralpha1model.m but uses a non-separable phase component.
%
%   beta            : STRF parameter vector
%                     beta(1): Response latency (msec)
%                     beta(2): Best Temporal Modulation Frequency (Hz)
%                     beta(3): Temporal modulation bandwidth (Hz)
%                     beta(4): Best octave frequency, xo
%                     beta(5): Best spectral modulation frequency (cycles/octaves)
%                     beta(6): Spectral Modulation Bandwidth (cycles/octaves) (sqrt(2)/(2*sigma))
%                     beta(7): Spectro-temporal phase (0-2*pi)
%                     beta(8): Peak envelope amplitude (excludes carrier)
%   input.taxis     : Time axis (msec)
%   input.X         : Octave frequency axis (octaves)
%
%RETURNED VARIABLES
%
%   STRFm           : Non-separable STRF model
%
% (C) F. He, M.A. Escabi, Feb/2018 (Edit 03/21; 08/22)
%
function [STRFm]=strfgaboralpha1modelns(beta,input)

    %Obtaining temporal (ET) and spectral (ES) envelopes
    %Envelops do not use phase information, which are only related to carriers
    betat=[beta(1:3) 0 1];
    betas=[beta(4:6) 0 1];  
    betas(3)=2*sqrt(2*log(2))/pi/betas(3);    % Convert spectral modulation bandwidth (cycles/oct) to Receptive Field Bandwidth (oct); Edit MAE/FH 08/22
    [~,ET]=trfalpha1model(betat,input.taxis);
    [~,ES]=srfgabormodel(betas,input.X);
    
    %Time and Frequency Axis
    taxis=input.taxis;
    X=input.X;
    X=(repmat(X',1,length(taxis)));
    taxis=repmat(taxis,size(X,1),1)/1000;
    
    %Computing nonseparable model STRF
    X0 = beta(4);
    RD = beta(5);
    Fm = beta(2);
    Del = beta(1)/1000;
    P = beta(7);
    STRFm = beta(8)*ET.*ES'.*cos(2*pi*RD*(X-X0)+2*pi*Fm*(taxis-Del)+P);     %The amplitude here is the envelope peak amplitude (not STRF peak amplitude because of carriers effects)

end
