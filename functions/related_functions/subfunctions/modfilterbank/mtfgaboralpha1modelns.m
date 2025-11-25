%function [H]=mtfgaboralpha1modelns(beta,input)
%
%   FILE NAME       : MTF GABOR ALPHA1 MODEL NS
%   DESCRIPTION     : This routine generates a non-separable MTF model where the MTF is the Fourier Transform
%                     of a separable gabor STRF. The model follows the same conventions as strfgabor1modelns.m
%                     but returned MTF is a frequency domain version of the STRF.
%
%                     The spectral receptive field is modeled as a gabor function (gaussian
%                     multiplied with cosine) while the temporal receptive field is modeled as the
%                     product of a 1-segment alpha fucntion and a cosine. The envelope of
%                     this model is seaparable, however, the carrier component
%                     (phase) is nonseparable. This allows the receptive fields to be selective
%                     for + and - spectro-temporal modulations.
%
%                     Note that the STRF from this model differ from strfgaboralpha1.m and
%                     strfgaboralpha2.m since they are generated using 1-segmenet alpha
%                     functions as opposed to 2-segment alpha functions.alpha. In those instances,
%                     1 and 2 designate the number of separable STRFs used in the model.
%
%                     This model is similar to strfgaboralpha1model.m but uses a non-separable phase component.
%
%   beta            : STRF parameter vector
%                     beta(1): Response latency (msec)
%                     beta(2): Fm: Character Modulation Frequency (Hz) (3dB)
%                     beta(3): Temporal modulation bandwidth (Hz)
%                     beta(4): Best octave frequency, xo
%                     beta(5): Spectral Modulation Bandwidth (cycles/octaves) (sqrt(2)/(2*sigma))
%                     beta(6): Best spectral modulation frequency (cycles/octaves)
%                     beta(7): Spectro-temporal phase (0-2*pi)
%                     beta(8): Peak envelope amplitude (excludes carrier)
%   input.FmAxis    : Modulation frequency axis (Hz)
%   input.RDAxis    : Ripple Density axis (cyc/oct)
%
%RETURNED VARIABLES
%
%   MTFm           : Non-separable MTF model
%
% (C) F. He, M.A. Escabi, Feb/2018 (Last Edit 03/21)
%
function [H]=mtfgaboralpha1modelns(beta,input)

    RDAxis = input.RDAxis;
    FmAxis = input.FmAxis;
    
%     delay = beta(1);    % onle delay = 0 is tested
%     x0 = beta(4);       % only x0 = 0 is tested
%     
%     Fm = beta(2);
%     bw_Fm = beta(3);
%     RD = beta(6);
    
    P = beta(7);
    bw_RD = beta(5);
    bw = 1/bw_RD;
    tau = 2*sqrt(sqrt(2)-1)/2/pi/beta(3);
    
    %Obtaining temporal (FmAxis) and spectral (RDAxis) envelopes
    phase_positive = 2*pi^2*exp(P*1i);
    SRF_freq_positive = bw/2*sqrt(pi)*exp(-bw^2/16*(2*pi*(RDAxis-beta(6))).^2);
    TRF_freq_positive = exp(1)*tau./(1+2*pi*tau*(FmAxis-beta(2))*1i).^2;
    phase_negative = 2*pi^2*exp(-P*1i);
    SRF_freq_negative = bw/2*sqrt(pi)*exp(-bw^2/16*(2*pi*(RDAxis+beta(6))).^2);
    TRF_freq_negative = exp(1)*tau./(1+2*pi*tau*(FmAxis+beta(2))*1i).^2;
    H = phase_positive * SRF_freq_positive' * TRF_freq_positive + phase_negative * SRF_freq_negative' * TRF_freq_negative;
    H = abs(H);
    H = H/max(max(H));

end

