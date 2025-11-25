%
%function [ModBankParam] = ModFiltBankParamGaborAlpha1ns(beta)
%
%   FILE NAME       : Modulation Transfer Function Filterbank
%                     GaborAlpha1 Non-separable
%
%   DESCRIPTION     : Generate arameters used to implement a STRF modulation
%                     filterbank using StrfModFiltBankGaborAlpha1ns.m. 
%                     
%                     The STRFs have a Gabor specral envelope and a 1st-order Alpha function
%                     temporal envelope and have a non-separable spectro-temporal phase. 
%                     
%   beta            : Parameter vector
%                     beta(1): Qt quality factor for temporal modualtion filters   (recommend: 1)
%                     beta(2): Qs quality factor for spectral modulation filters  (recommend: 1)
%                     beta(3): Fml - lower temporal modulation frequency (Hz)
%                     beta(4): Fmu - upper temporal modulation frequency (Hz)
%                     beta(5): RDl - lower spectral modulation frequency (cycles/oct)
%                     beta(6): RDu - upper spectral modulation frequency (cycles/oct)
%                     beta(7): Dt  - temporal modulation filter spacing (octave)
%                     beta(8): Ds  - spectral modulation filter spacing (octave)
%
%RETURNED VARIABLES
%
%   ModBankParam.F(j,i) - Matrix of data structure containg Modulation Filterbank STRF parameters 
%
%                .Delay : Response latency (msec)
%                .Fm    : Character Modulation Frequency (Hz) (3dB)
%                .BWt   : Temporal modulation bandwidth (Hz)
%                .Bof   : Best octave frequency, xo (default = 0)
%                .RD    : Best spectral modulation frequency (cycles/octaves)
%                .BWs   : Spectral Modulation Bandwidth (cycles/octaves) (sqrt(2)/(2*sigma))
%                .P     : Spectro-temporal phase (0-2*pi) (default pi/4)
%                .Amp   : Amplitude (excludes carrier)
%                .Beta  : Everything above in sequential order
%   
%   ModBankParam.Param - Modulation Filter Bank Parameter Data Structure
%
%                .Qt     : Qt quality factor for temporal modualtion filters
%                .Qs     : Qs quality factor for spectral modualtion filters
%                .Fml    : Lower temporal modulation frequency (Hz)
%                .Fmu    : Upper temporal modulation frequency (Hz)
%                .RDl    : Lower spectral modulation frequency (cycles/oct
%                .RDu    : Upper spectral modulation frequency (cycles/oct)
%                .Dt     : Temporal modulation filter spacing (octave)
%                .Ds     : Spectral modulation filter spacing (octave)
%                .beta   : Original input vector containing filterbank parameters
%
% C) F. He & M.A. Escabi, April 2021 (Edit Oct 2021, MAE)
%
function [ModBankParam] = ModFiltBankParamGaborAlpha1ns(beta)

%Modulation Filter Parameters
Qt = beta(1);
Qs = beta(2);
Fml = beta(3);
Fmu = beta(4);
RDl = beta(5);
RDu = beta(6);
Dt = beta(7);   %Temporal modulatin filter spacing in octaves
Ds = beta(8);   %Spectral modulatin filter spacing in octaves

%Generate character modulation frequency sequence
%temporal
TN=log2(Fmu/Fml);
L=ceil(TN/Dt);
Tc=(0:L)*Dt;
Fm=Fml*2.^Tc;
Fm = [0 Fm];

%spectral
SN=log2(RDu/RDl);
L=ceil(SN/Ds);
Sc=(0:L)*Ds;
RD=RDl*2.^Sc;
RD = [0 RD];

%Generate modulation filterbank parameters
Delay = 0;                      %value other than 0 will cause problems
Tb = Fm./Qt;                    %temporal bandwidth (Hz)
Tb(1) = 2*Fm(2);                %Bandwidth of DC filter - note that it is 2 x the Fm of the first non DC filter
% Tphase = pi/4;                  %Temporal phase

Bof = 0;                        %Best Octave Frequency
Sb = RD./Qs;                    %spectral bandwidth (cycyles/octive Hz)
Sb(1) = 2*RD(2);                %Bandwidth of DC filter - note that it is 2 x the RD of the first non DC filter
% Sphase = 0;                     %Spectral phase
Amp = 1;                        %Amplitude

%Generate negative Fm filters
Fm = [-flip(Fm(2:end)) Fm];
Tb = [flip(Tb(2:end)) Tb];
P = -pi/4;                       %Spectro-temporal phase, consistent with the seperable versions

%Generate Modualtion filters
Nt = length(Fm);                %Number of temporal filters
Ns = length(RD);                %Number of spectral filters
for i = 1:Nt
    for j = 1:Ns
        if Fm(i)>0
            Beta=[Delay Fm(i) Tb(i) Bof RD(j) Sb(j) P Amp];
        else
            if Fm(i) == 0
                Beta=[Delay Fm(i) Tb(i) Bof RD(j) Sb(j) 0 Amp];
            else
                Beta=[Delay Fm(i) Tb(i) Bof RD(j) Sb(j) -P Amp];
            end
        end
        
        %STRF Filter parameters
        ModBankParam.F(j,i).Delay = Beta(1);
        ModBankParam.F(j,i).Fm = Beta(2);
        ModBankParam.F(j,i).BWt = Beta(3);
        ModBankParam.F(j,i).Bof = Beta(4);
        ModBankParam.F(j,i).RD = Beta(5);
        ModBankParam.F(j,i).BWs = Beta(6);
        ModBankParam.F(j,i).P = Beta(7);
        ModBankParam.F(j,i).Amp = Beta(8);
        ModBankParam.F(j,i).Beta = Beta;                            %Contains the Gabor Alpha1 STRF Parameters
    end
end

%Modualtion Filterbank Paramneters
ModBankParam.Param.Qt        = beta(1);
ModBankParam.Param.Qs        = beta(2);
ModBankParam.Param.Fml       = beta(3);
ModBankParam.Param.Fmu       = beta(4);
ModBankParam.Param.RDl       = beta(5);
ModBankParam.Param.RDu       = beta(6);
ModBankParam.Param.Dt        = beta(7);
ModBankParam.Param.Ds        = beta(8);
ModBankParam.Param.beta      = beta;                                 %Original Beta Parameter Vector (Filterbank Parameters)
ModBankParam.Param.FmAxis    = [ModBankParam.F(1,:).Fm];
ModBankParam.Param.RDAxis    = [ModBankParam.F(:,1).RD];
ModBankParam.Param.LogFmAxis = [-fliplr(Tc),-Inf,Tc];                %Give back the log2 scale FmAxis (with Fm(2) as reference)
ModBankParam.Param.LogRDAxis = [-Inf,Sc];
end
