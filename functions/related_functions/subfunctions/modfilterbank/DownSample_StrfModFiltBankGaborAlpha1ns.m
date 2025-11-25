%
%function [StrfBank] = StrfModFiltBankGaborAlpha1ns(beta,axis,Fses,cut_th)
%
%   FILE NAME       : Spectral Temporal Receptive Field Filterbank
%                     GaborAlpha1 Non-separable
%
%   DESCRIPTION     : Generate the whole filter bank and return H
%
%   beta            : Parameter vector
%                     beta(1): Qt quality factor for temporal   (recommend: 1)
%                     beta(2): Qs quality factor for spectral   (recommend: 1)
%                     beta(3): Fml temporal lower CF
%                     beta(4): Fmu temporal upper CF 
%                     beta(5): RDl spectral lower CF (ripple density)
%                     beta(6): RDu spectral upper CF (ripple density)
%                     beta(7): Dt temporal resolution (octave)
%                     beta(8): Ds spectral resolution (octave)
%
%   axis            : axis.taxis in msec
%                     axis.X in cyc/oct
%   
%   Fses            : sample frequencies for axis
%                     Fses(1) for temporal
%                     Fses(2) for spectral
%                     Note: Fses and axis are exclusive (they describe the same thing in different ways)
% 
%   cut_th          : cut off threshold for filters (default 0.01, one percent)
% 
%RETURNED VARIABLES
%
%   StrfBank.F      : StrfBank.F(j,i).Delay : Response latency (msec)
%                     StrfBank.F(j,i).Fm : Character Modulation Frequency (Hz) (3dB)
%                     StrfBank.F(j,i).BWt : Temporal modulation bandwidth (Hz)
%                     StrfBank.F(j,i).Bof : Best octave frequency, xo (default = 0)
%                     StrfBank.F(j,i).BWs : Spectral Modulation Bandwidth (cycles/octaves) (sqrt(2)/(2*sigma))
%                     StrfBank.F(j,i).RD : Best spectral modulation frequency (cycles/octaves)
%                     StrfBank.F(j,i).P : Spectro-temporal phase (0-2*pi) (default pi/4)
%                     StrfBank.F(j,i).Amp : Amplitude (excludes carrier)
%                     StrfBank.F(j,i).Beta : Everything above
%                     StrfBank.F(j,i).input : axis for the filter
%                     StrfBank.F(j,i).H : filter
%
%   StrfBank.Param  : Keep the input beta
%
%NOTE
%   axis            : if not manually given, default is 2000msec for time, +/-4cyc/oct for spectral
%                     suitable for filterbank with lowest frequency ~ (1Hz, 0.25cyc/oct) or DC with similar bw
% FH Apr/2021

function [StrfBank] = StrfModFiltBankGaborAlpha1ns(beta,axis,Fses,cut_th)
if nargin<4 || isempty(cut_th)
    cut_th = 0.01;
end 
if nargin<3 || isempty(Fses)
    if nargin<2 || isempty(axis)
        Fsesflag = 0;
        axisflag = 0;
    else
        Fsesflag = 0;
        axisflag = 1;
    end
else
    if isempty(axis)
        Fsesflag = 1;
        axisflag = 0;
    else
        Fsesflag = 1;
        axisflag = 1;
    end
end

%Modulation Filter Parameters
Qt = beta(1);
Qs = beta(2);
Fml = beta(3);
Fmu = beta(4);
RDl = beta(5);
RDu = beta(6);
Dt = beta(7);
Ds = beta(8);

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
P = -pi/4;                       %Spectro-temporal phase !!! inorder to be consistent with seperable version

%Generate axis(input)
TMAX = 2; % second
XMAX = 4; 
if Fsesflag
    if axisflag
        warning('axis has two relative informations, which may conflict');
    else
        dt = 1/Fses(1);
        tmax = TMAX;
        dX = 1/Fses(2);
        Xmax = XMAX;
        input.taxis = (0:dt:tmax)*1000;
        input.X = -Xmax:dX:Xmax;
    end
else
    if axisflag
        input = axis;
    else
        dt = 1/2/(max(Fm)+max(Tb)); % highest frequency = max(Fm)+1/2*max(Tb)
        tmax = TMAX;
        dX = 1/2/(max(RD)+max(Sb)); % highest ripple density = max(RD)+1/2*max(Sb)
        Xmax = XMAX;
        input.taxis = (0:dt:tmax)*1000; % msec
        input.X = -Xmax:dX:Xmax;
    end
end

%Generate Modualtion filters
Nt = length(Fm);                %Number of temporal filters
Ns = length(RD);                %Number of spectral filters

%% H
% tDC and sDC
i = (Nt-1)/2+1;
j = 1;
Beta=[Delay Fm(i) Tb(i) Bof Sb(j) RD(j) 0 Amp];
H = strfgaboralpha1modelns(Beta,input);
% truncate filter and axis
axisx = [0,0];
axisy = [0,0];
temp_H = abs(sum(H,1));
axisx(1) = 1;
axisx(2) = find(temp_H>cut_th*max(temp_H),1,'last');
temp_H = abs(sum(H,2));
axisy(1) = find(temp_H>cut_th*max(temp_H),1,'first');
axisy(2) = find(temp_H>cut_th*max(temp_H),1,'last');
temp_input.taxis = input.taxis(axisx(1):axisx(2));
temp_input.X = input.X(axisy(1):axisy(2));
H = H(axisy(1):axisy(2),axisx(1):axisx(2));
% save to output
StrfBank.F(j,i).Delay = Delay;
StrfBank.F(j,i).Fm = Fm(i);
StrfBank.F(j,i).BWt = Tb(i);
StrfBank.F(j,i).Bof = Bof;
StrfBank.F(j,i).BWs = Sb(j);
StrfBank.F(j,i).RD = RD(j);
StrfBank.F(j,i).P = P;
StrfBank.F(j,i).Amp = Amp;
StrfBank.F(j,i).Beta = Beta;
StrfBank.F(j,i).input = temp_input;
StrfBank.F(j,i).H = H;

% sDC
j = 1;
i = (Nt-1)/2+2;
Beta = [Delay Fm(i) Tb(i) Bof Sb(j) RD(j) P Amp];
H = strfgaboralpha1modelns(Beta,input);
% truncate filter and axis
axisx = [0,0];
axisy = [0,0];
temp_H = abs(sum(H,1));
axisx(1) = 1;
axisx(2) = find(temp_H>cut_th*max(temp_H),1,'last');
temp_H = abs(sum(H,2));
axisy(1) = find(temp_H>cut_th*max(temp_H),1,'first');
axisy(2) = find(temp_H>cut_th*max(temp_H),1,'last');
temp_input.taxis = input.taxis(axisx(1):axisx(2));
temp_input.X = input.X(axisy(1):axisy(2));
H = H(axisy(1):axisy(2),axisx(1):axisx(2));
% save to output
StrfBank.F(j,i).Delay = Delay;
StrfBank.F(j,i).Fm = Fm(i);
StrfBank.F(j,i).BWt = Tb(i);
StrfBank.F(j,i).Bof = Bof;
StrfBank.F(j,i).BWs = Sb(j);
StrfBank.F(j,i).RD = RD(j);
StrfBank.F(j,i).P = P;
StrfBank.F(j,i).Amp = Amp;
StrfBank.F(j,i).Beta = Beta;
StrfBank.F(j,i).input = temp_input;
StrfBank.F(j,i).H = H;
% opposite filter save to output
StrfBank.F(j,i-2).Delay = Delay;
StrfBank.F(j,i-2).Fm = Fm(i);
StrfBank.F(j,i-2).BWt = Tb(i);
StrfBank.F(j,i-2).Bof = Bof;
StrfBank.F(j,i-2).BWs = Sb(j);
StrfBank.F(j,i-2).RD = RD(j);
StrfBank.F(j,i-2).P = -P;
StrfBank.F(j,i-2).Amp = Amp;
StrfBank.F(j,i-2).Beta = Beta;
StrfBank.F(j,i-2).input = temp_input;
StrfBank.F(j,i-2).H = H;
for i = (Nt-1)/2+3:Nt
    Beta = [Delay Fm(i) Tb(i) Bof Sb(j) RD(j) P Amp];
    % save to output
    StrfBank.F(j,i).Delay = Delay;
    StrfBank.F(j,i).Fm = Fm(i);
    StrfBank.F(j,i).BWt = Tb(i);
    StrfBank.F(j,i).Bof = Bof;
    StrfBank.F(j,i).BWs = Sb(j);
    StrfBank.F(j,i).RD = RD(j);
    StrfBank.F(j,i).P = P;
    StrfBank.F(j,i).Amp = Amp;
    StrfBank.F(j,i).Beta = Beta;
    StrfBank.F(j,i).H = StrfBank.F(j,i-1).H(:,1:2:end);
    StrfBank.F(j,i).input.taxis = StrfBank.F(j,i-1).input.taxis(1:size(StrfBank.F(j,i).H,2));
    StrfBank.F(j,i).input.X = StrfBank.F(j,i-1).input.X;
end
for i = (Nt-1)/2-1:-1:1
    Beta = [Delay Fm(i) Tb(i) Bof Sb(j) RD(j) -P Amp];
    % save to output
    StrfBank.F(j,i).Delay = Delay;
    StrfBank.F(j,i).Fm = Fm(i);
    StrfBank.F(j,i).BWt = Tb(i);
    StrfBank.F(j,i).Bof = Bof;
    StrfBank.F(j,i).BWs = Sb(j);
    StrfBank.F(j,i).RD = RD(j);
    StrfBank.F(j,i).P = -P;
    StrfBank.F(j,i).Amp = Amp;
    StrfBank.F(j,i).Beta = Beta;
    StrfBank.F(j,i).H = StrfBank.F(j,i+1).H(:,1:2:end);
    StrfBank.F(j,i).input.taxis = StrfBank.F(j,i+1).input.taxis(1:size(StrfBank.F(j,i).H,2));
    StrfBank.F(j,i).input.X = StrfBank.F(j,i+1).input.X;
end

% tDC
j = 2;
i = (Nt-1)/2+1;
Beta = [Delay Fm(i) Tb(i) Bof Sb(j) RD(j) 0 Amp];
H = strfgaboralpha1modelns(Beta,input);
% truncate filter and axis
axisx = [0,0];
axisy = [0,0];
temp_H = abs(sum(H,1));
axisx(1) = 1;
axisx(2) = find(temp_H>cut_th*max(temp_H),1,'last');
temp_H = abs(sum(H,2));
axisy(1) = find(temp_H>cut_th*max(temp_H),1,'first');
axisy(2) = find(temp_H>cut_th*max(temp_H),1,'last');
temp_input.taxis = input.taxis(axisx(1):axisx(2));
temp_input.X = input.X(axisy(1):axisy(2));
H = H(axisy(1):axisy(2),axisx(1):axisx(2));
% save to output
StrfBank.F(j,i).Delay = Delay;
StrfBank.F(j,i).Fm = Fm(i);
StrfBank.F(j,i).BWt = Tb(i);
StrfBank.F(j,i).Bof = Bof;
StrfBank.F(j,i).BWs = Sb(j);
StrfBank.F(j,i).RD = RD(j);
StrfBank.F(j,i).P = 0;
StrfBank.F(j,i).Amp = Amp;
StrfBank.F(j,i).Beta = Beta;
StrfBank.F(j,i).input = temp_input;
StrfBank.F(j,i).H = H;
for j = 3:Ns
    Beta = [Delay Fm(i) Tb(i) Bof Sb(j) RD(j) 0 Amp];
    lenX = length(StrfBank.F(j-1,i).input.X);
    % save to output
    StrfBank.F(j,i).Delay = Delay;
    StrfBank.F(j,i).Fm = Fm(i);
    StrfBank.F(j,i).BWt = Tb(i);
    StrfBank.F(j,i).Bof = Bof;
    StrfBank.F(j,i).BWs = Sb(j);
    StrfBank.F(j,i).RD = RD(j);
    StrfBank.F(j,i).P = 0;
    StrfBank.F(j,i).Amp = Amp;
    StrfBank.F(j,i).Beta = Beta;
    StrfBank.F(j,i).H = [flip(StrfBank.F(j-1,i).H((lenX-1)/2+1:-2:1,:),1);StrfBank.F(j-1,i).H((lenX-1)/2+1:2:end,:)];
    StrfBank.F(j,i).input.taxis = StrfBank.F(j-1,i).input.taxis;
    StrfBank.F(j,i).input.X = StrfBank.F(j-1,i).input.X((lenX-1)/2+1-floor((lenX-1)/4):(lenX-1)/2+1+floor((lenX-1)/4));
end

% others
j = 2;
i = (Nt-1)/2+2;
Beta = [Delay Fm(i) Tb(i) Bof Sb(j) RD(j) P Amp];
H = strfgaboralpha1modelns(Beta,input);
% truncate filter and axis
axisx = [0,0];
axisy = [0,0];
temp_H = abs(sum(H,1));
axisx(1) = 1;
axisx(2) = find(temp_H>cut_th*max(temp_H),1,'last');
temp_H = abs(sum(H,2));
axisy(1) = find(temp_H>cut_th*max(temp_H),1,'first');
axisy(2) = find(temp_H>cut_th*max(temp_H),1,'last');
temp_input.taxis = input.taxis(axisx(1):axisx(2));
temp_input.X = input.X(axisy(1):axisy(2));
H = H(axisy(1):axisy(2),axisx(1):axisx(2));
% save to output
StrfBank.F(j,i).Delay = Delay;
StrfBank.F(j,i).Fm = Fm(i);
StrfBank.F(j,i).BWt = Tb(i);
StrfBank.F(j,i).Bof = Bof;
StrfBank.F(j,i).BWs = Sb(j);
StrfBank.F(j,i).RD = RD(j);
StrfBank.F(j,i).P = P;
StrfBank.F(j,i).Amp = Amp;
StrfBank.F(j,i).Beta = Beta;
StrfBank.F(j,i).input = temp_input;
StrfBank.F(j,i).H = H;
% opposite filter save to output
i = (Nt-1)/2;
Beta = [Delay Fm(i) Tb(i) Bof Sb(j) RD(j) -P Amp];
H = strfgaboralpha1modelns(Beta,input);
% truncate filter and axis
axisx = [0,0];
axisy = [0,0];
temp_H = abs(sum(H,1));
axisx(1) = 1;
axisx(2) = find(temp_H>cut_th*max(temp_H),1,'last');
temp_H = abs(sum(H,2));
axisy(1) = find(temp_H>cut_th*max(temp_H),1,'first');
axisy(2) = find(temp_H>cut_th*max(temp_H),1,'last');
temp_input.taxis = input.taxis(axisx(1):axisx(2));
temp_input.X = input.X(axisy(1):axisy(2));
H = H(axisy(1):axisy(2),axisx(1):axisx(2));
% save to output
StrfBank.F(j,i).Delay = Delay;
StrfBank.F(j,i).Fm = Fm(i);
StrfBank.F(j,i).BWt = Tb(i);
StrfBank.F(j,i).Bof = Bof;
StrfBank.F(j,i).BWs = Sb(j);
StrfBank.F(j,i).RD = RD(j);
StrfBank.F(j,i).P = -P;
StrfBank.F(j,i).Amp = Amp;
StrfBank.F(j,i).Beta = Beta;
StrfBank.F(j,i).input = temp_input;
StrfBank.F(j,i).H = H;
for j = 3:Ns
    % save to output
    i = (Nt-1)/2+2;
    lenX = length(StrfBank.F(j-1,i).input.X);
    Beta = [Delay Fm(i) Tb(i) Bof Sb(j) RD(j) P Amp];
    StrfBank.F(j,i).Delay = Delay;
    StrfBank.F(j,i).Fm = Fm(i);
    StrfBank.F(j,i).BWt = Tb(i);
    StrfBank.F(j,i).Bof = Bof;
    StrfBank.F(j,i).BWs = Sb(j);
    StrfBank.F(j,i).RD = RD(j);
    StrfBank.F(j,i).P = P;
    StrfBank.F(j,i).Amp = Amp;
    StrfBank.F(j,i).Beta = Beta;
    StrfBank.F(j,i).H = [flip(StrfBank.F(j-1,i).H((lenX-1)/2+1:-2:1,:),1);StrfBank.F(j-1,i).H((lenX-1)/2+1:2:end,:)];
    StrfBank.F(j,i).input.taxis = StrfBank.F(j-1,i).input.taxis;
    StrfBank.F(j,i).input.X = StrfBank.F(j-1,i).input.X((lenX-1)/2+1-floor((lenX-1)/4):(lenX-1)/2+1+floor((lenX-1)/4));
    % opposite filter save to output
    i = (Nt-1)/2;
    lenX = length(StrfBank.F(j-1,i).input.X);
    Beta = [Delay Fm(i) Tb(i) Bof Sb(j) RD(j) -P Amp];
    StrfBank.F(j,i).Delay = Delay;
    StrfBank.F(j,i).Fm = Fm(i);
    StrfBank.F(j,i).BWt = Tb(i);
    StrfBank.F(j,i).Bof = Bof;
    StrfBank.F(j,i).BWs = Sb(j);
    StrfBank.F(j,i).RD = RD(j);
    StrfBank.F(j,i).P = -P;
    StrfBank.F(j,i).Amp = Amp;
    StrfBank.F(j,i).Beta = Beta;
    StrfBank.F(j,i).H = [flip(StrfBank.F(j-1,i).H((lenX-1)/2+1:-2:1,:),1);StrfBank.F(j-1,i).H((lenX-1)/2+1:2:end,:)];
    StrfBank.F(j,i).input.taxis = StrfBank.F(j-1,i).input.taxis;
    StrfBank.F(j,i).input.X = StrfBank.F(j-1,i).input.X((lenX-1)/2+1-floor((lenX-1)/4):(lenX-1)/2+1+floor((lenX-1)/4));
end
for j = 3:Ns
    for i = (Nt-1)/2+3:Nt
        Beta = [Delay Fm(i) Tb(i) Bof Sb(j) RD(j) P Amp];
        % save to output
        StrfBank.F(j,i).Delay = Delay;
        StrfBank.F(j,i).Fm = Fm(i);
        StrfBank.F(j,i).BWt = Tb(i);
        StrfBank.F(j,i).Bof = Bof;
        StrfBank.F(j,i).BWs = Sb(j);
        StrfBank.F(j,i).RD = RD(j);
        StrfBank.F(j,i).P = P;
        StrfBank.F(j,i).Amp = Amp;
        StrfBank.F(j,i).Beta = Beta;
        StrfBank.F(j,i).H = StrfBank.F(j,i-1).H(:,1:2:end);
        StrfBank.F(j,i).input.taxis = StrfBank.F(j,i-1).input.taxis(1:size(StrfBank.F(j,i).H,2));
        StrfBank.F(j,i).input.X = StrfBank.F(j,i-1).input.X;
    end
end
for j = 3:Ns
    for i = (Nt-1)/2-1:-1:1
        Beta = [Delay Fm(i) Tb(i) Bof Sb(j) RD(j) -P Amp];
        % save to output
        StrfBank.F(j,i).Delay = Delay;
        StrfBank.F(j,i).Fm = Fm(i);
        StrfBank.F(j,i).BWt = Tb(i);
        StrfBank.F(j,i).Bof = Bof;
        StrfBank.F(j,i).BWs = Sb(j);
        StrfBank.F(j,i).RD = RD(j);
        StrfBank.F(j,i).P = -P;
        StrfBank.F(j,i).Amp = Amp;
        StrfBank.F(j,i).Beta = Beta;
        StrfBank.F(j,i).H = StrfBank.F(j,i+1).H(:,1:2:end);
        StrfBank.F(j,i).input.taxis = StrfBank.F(j,i+1).input.taxis(1:size(StrfBank.F(j,i).H,2));
        StrfBank.F(j,i).input.X = StrfBank.F(j,i+1).input.X;
    end
end

%% Modualtion Filterbank Paramneters
StrfBank.Param.Qt        = beta(1);
StrfBank.Param.Qs        = beta(2);
StrfBank.Param.Fml       = beta(3);
StrfBank.Param.Fmu       = beta(4);
StrfBank.Param.RDl       = beta(5);
StrfBank.Param.RDu       = beta(6);
StrfBank.Param.Dt        = beta(7);
StrfBank.Param.Ds        = beta(8);
StrfBank.Param.beta      = beta;                                 %Original Beta Vector
StrfBank.Param.cut_th    = cut_th;
StrfBank.Param.input     = input;
end
