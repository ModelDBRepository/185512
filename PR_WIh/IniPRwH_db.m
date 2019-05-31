function aPRwH=IniPRwH_db(Isinj,Ek,gh,h_Vhalf,varargin)

%% May 23 2015
%% June 11 2015 added Units, Rdsin, Memb Area
% First load in defaults based on Pinsky-Rinzel 1994
% and Isinj and Ek
% Then overwrite paramters depending on varargin
        aPRwH.Cm = 3; % 10^-6 F/cm^2 (uF/cm^2)
        aPRwH.gL=0.1; %10^-3 S/cm^2 (mS/cm^2)
        aPRwH.gNa=30.0; %mS/cm^2
        aPRwH.gKDR=15.0; %mS/cm^2
        aPRwH.gKC = 15.0; %mS/cm^2
        aPRwH.gKAHP = 0.8; %mS/cm^2
        aPRwH.gCa = 10.0; % mS/cm^2
        aPRwH.gh= gh; %0.03; %mS/cm^2 %control level Lippert 2009
        aPRwH.h_Vhalf=h_Vhalf;% -81; % mV control level Lippert 2009
        aPRwH.ENa = 120.0; %mV WRT - 60 mV (80 real or Booth and Bose mV)
        aPRwH.Ek=Ek; 
        aPRwH.EL = 0; % mV
        aPRwH.ECa = 140.0; % mV
        aPRwH.Eh = 35.0; %-25 + 60; % mV Lippert 2009 kept it constnat Dyrhfjeld_Johnson 2009 refernces Robinson and Siegelbaum, 2003 a range of -25 mV to -40 mV) Not sure why Lippert picked the high end?
        aPRwH.p = 0.5; %unitless fraction area soma
        aPRwH.gc = 2.1; %mS/cm^2
        aPRwH.WRT = -60;
        aPRwH.HcurrLippertWRT=0;
        aPRwH.DiffRefVoltHcurr=aPRwH.HcurrLippertWRT-aPRwH.WRT;
        aPRwH.Vsyn=60;
        aPRwH.MaxS=125;
        aPRwH.Isinj=Isinj;
        aPRwH.Idinj=0;
        aPRwH.TotMembArea = 6e-6; %cm^2
        aPRwH.Units.AreaScaleToM = 1e-4;
        aPRwH.Units.VoltsScaleToV=1e-3;
        aPRwH.Units.CondPerUnitAreaToSpercmsqr=1e-3;
        aPRwH.Units.InjCurrentToAmpsPercmsqr=1e-6;
        aPRwH.Units.Notes={'all conductances including NMDA and AMPA are in mS/cm^2','All voltages are in mV'};
        aPRwH.RDSinOhms=1/(aPRwH.gc*aPRwH.Units.CondPerUnitAreaToSpercmsqr*aPRwH.TotMembArea);
        aPRwH.RDSinMOhms=aPRwH.RDSinOhms/1e+6;
  %      aPRSS = GetSinglePRSSTonic(1,1000,aPRwH.Idinj,aPRwH.Isinj);
  %      aPRwH.QuasiSS=aPRSS;
        aPRwH.Units.AreaScaleToM = 1e-4;
        aPRwH.Units.VoltsScaleToV=1e-3;
        aPRwH.Units.CondPerUnitAreaToSpercmsqr=1e-3;
        aPRwH.Units.InjCurrentToAmpsPercmsqr=1e-6;
        aPRwH.Units.Notes={'all conductances including NMDA and AMPA are in mS/cm^2','All voltages are in mV'};

if nargin > 9
    display('too many variable arguments')
    exit
end
if nargin>=5
    aPRwH.ENa=varargin{1};
end
if nargin>=6
    aPRwH.p=varargin{2};
end
if nargin>=75
    aPRwH.gc=varargin{3};
end
if nargin>=8
    aPRwH.gKAHP=varargin{4};
end
if nargin>=7
    aPRwH.gKC=varargin{5};
if nargin>=8
    aPRwH.Idinj=varargin{6};
end

end