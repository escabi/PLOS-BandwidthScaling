%function [ModFiltBank] = ModFiltBankBSpline(beta,Fsd,hflag)
%
%       FILE NAME       : MOD FILTBANK BSPLINE
%       DESCRIPTION     : Generate a B-Spline modualtion filterbank
%
%       beta            : Parameter vector
%                         beta(1): Qt quality factor for temporal
%                         beta(2): Qs quality factor for spectral
%                         beta(3): Fml temporal lower Mod Frequency
%                         beta(4): Fmu temporal upper Mod Frequency
%                         beta(5): RDl spectral lower Mod Frequency (ripple density)
%                         beta(6): RDu spectral upper Mod Frequency (ripple density)
%       Fsd             : Vector containing the desired sampling rate for
%                         the temporal and spectral axis
%
%                           Fsd = [Fsd_t Fsd_s]
%
%                         where Fsd_t is in Hz and Fsd_s is in
%                         samples/octave.
%       hflag           : 1 or 2 - Optional (Default == 2)
%                         1 - Return Gabor IR using optimal IR matrix size for
%                             each scale
%                         2 - Return Gabor IR using a standard IR matrix size
%
%RETURNED VARIABLES
%
%       ModFildBank     : Data strcture containg modulation filterbank 
%                         parameters and filters
%
%           .F          : Vector of data structures containg impulse
%                         response of the B-spline filters
%           .Param      : Data structure cotaining the filterbank
%                         parameters
%
% (C) MAE, 03/20
%
function [ModFiltBank] = ModFiltBankBSpline(beta,Fsd,hflag)

%Input Args
if nargin<3 || isempty(hflag)
    hflag = 3;
end

%Modulation Filter Parameters
Qt = beta(1);
Qs = beta(2);
Fml = beta(3);
Fmu = beta(4);
RDl = beta(5);
RDu = beta(6);

%Sampling Rates
Fsd_t=Fsd(1);
Fsd_s=Fsd(2);

%Temporal modulation frequency sequence
TN=log2(Fmu/Fml);                   %Number of octaves
Dt=log2((2*Qt+1)/(2*Qt-1));         %Filter resolution in octaves
L=ceil(TN/Dt);                      %Number of filters - excluding DC
Fml=Fmu/((2*Qt+1)/(2*Qt-1))^L;      %Recomputing lower temporal mod frequency - make sure upper is exact
Fm=Fml*((2*Qt+1)/(2*Qt-1)).^(0:L);
Fm=[0 Fm];                          %Temporal Modualtion Cutoff Sequence
BWt=Fm(2:end)-Fm(1:end-1);          %Temporal Modulation Bandwidth
BWt(1)=2*Fm(2);                     %Temporal Modualtion Bandwidth of DC filter - includes negative frequencies

%Spectral modulation frequency sequence
SN=log2(RDu/RDl);                   %Number of octaves
Ds=log2((2*Qs+1)/(2*Qs-1));         %Filter resolution in octaves
L=ceil(SN/Ds);                      %Number of filters - excluding DC
RDl=RDu/((2*Qs+1)/(2*Qs-1))^L;      %Recomputing lower spectral mod frequency   make sure upper is exact
RD=RDl*((2*Qs+1)/(2*Qs-1)).^(0:L);
RD=[0 RD];                          %Spectral Modualtion Cutoff Sequence
BWs=RD(2:end)-RD(1:end-1);          %Spectral Modulation Bandwidth
BWs(1)=2*RD(2);                     %Spectral Modualtion Bandwidth of DC filter - includes negative frequencies

%Temporal B-Spline Modualtion Filterbank
ATT=30;
for k=1:length(Fm)-1
    if Fm(k)==0
        [Ht] = lowpass(BWt(k)/2,BWt(k)/2,Fsd_t,ATT);   
    else
        [Ht] = bandpass(Fm(k),Fm(k+1),BWt(k)/2,Fsd_t,ATT);    
    end
    Nt=(length(Ht)-1)/2;
    if k==1
        i=min(find(min(BWt)==BWt));
        [HH] = bandpass(Fm(i),Fm(i+1),BWt(i)/2,Fsd_t,ATT);
        N=(length(HH)-1)/2;
    end
    
    %Selecting FIR filter size
    if hflag==2
        FiltTemp(k).H=zeros(size(HH));
        FiltTemp(k).H(N+1+(-Nt:Nt))=Ht;
    elseif hflag==1
        FiltTemp(k).H=Ht;
    end
end

%Spectral B-Spline Modulation Filterbank
ATT=30;
for k=1:length(RD)-1
    if Fm(k)==0
        [Hs] = lowpass(BWs(k)/2,BWs(k)/2,Fsd_s,ATT);   
    else
        [Hs] = bandpass(RD(k),RD(k+1),BWs(k)/2,Fsd_s,ATT);    
    end
    Ns=(length(Hs)-1)/2;
    if k==1
        i=min(find(min(BWs)==BWs));
        [HH] = bandpass(RD(i),RD(i+1),BWs(i)/2,Fsd_s,ATT);  
        N=(length(HH)-1)/2;
    end
    
    %Selecting FIR filter size
    if hflag==2
        FiltSpec(k).H=zeros(size(HH));
        FiltSpec(k).H(N+1+(-Ns:Ns))=Hs;
    elseif hflag==1
        FiltSpec(k).H=Hs;
    end
end

%Generating Spectro-temporal Impulse Responses
for k=1:length(RD)-1
    for l=1:length(Fm)-1
        ModFiltBank.F(k,l).H=FiltSpec(k).H'*FiltTemp(l).H;
    end
end

% 
% For testing
% H=zeros(size(ModFiltBank.F(1,1).H));
% for k=1:size(ModFiltBank.F,1)-1
%     for l=1:size(ModFiltBank.F,2)-1
%         subplot(5,7,l+(k-1)*7)
%         imagesc(ModFiltBank.F(k,l).H)
%         H=H+ModFiltBank.F(k,l).H;
%     end
% end
% figure
% imagesc(20*log10(fftshift(abs(fft2(H,1024*8,1024*8))))),colorbar,colormap jet,caxis([-60 3])

%Modualtion Filterbank Paramneters
ModFiltBank.Param.Qt        = beta(1);
ModFiltBank.Param.Qs        = beta(2);
ModFiltBank.Param.Fml       = beta(3);
ModFiltBank.Param.Fmu       = beta(4);
ModFiltBank.Param.RDl       = beta(5);
ModFiltBank.Param.RDu       = beta(6);
ModFiltBank.Param.Dt        = Dt;
ModFiltBank.Param.Ds        = Ds;
ModFiltBank.Param.Fm        = Fm;
ModFiltBank.Param .BWt      = BWt;
ModFiltBank.Param.RD        = RD;
ModFiltBank.Param .BWs      = BWs;
ModFiltBank.Param.Fsd_t     = Fsd(1);
ModFiltBank.Param.Fsd_s     = Fsd(2);
ModFiltBank.Param.beta      = beta;                                 %Original Beta Vector