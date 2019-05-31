function aPR=IniPR_db(Isinj,Ek,varargin)

%% May 23 2015
%% June 11 2015 added Units, Rdsin, Memb Area
% First load in defaults based on Pinsky-Rinzel 1994
% and Isinj and Ek
% Then overwrite paramters depending on varargin
        aPR.Cm = 3; % 10^-6 F/cm^2 (uF/cm^2)
        aPR.gL=0.1; %10^-3 S/cm^2 (mS/cm^2)
        aPR.gNa=30.0; %mS/cm^2
        aPR.gKDR=15.0; %mS/cm^2
        aPR.gKC = 15.0; %mS/cm^2
        aPR.gKAHP = 0.8; %mS/cm^2
        aPR.gCa = 10.0; % mS/cm^2
        aPR.ENa = 120.0; %mV WRT - 60 mV (80 real or Booth and Bose mV)
        aPR.Ek=Ek; 
        aPR.EL = 0; % mV
        aPR.ECa = 140.0; % mV
        aPR.p = 0.5; %unitless fraction area soma
        aPR.gc = 2.1; %mS/cm^2
        aPR.WRT = -60;
        aPR.Vsyn=60;
        aPR.MaxS=125;
        aPR.Isinj=Isinj;
        aPR.Idinj=0;
        aPR.TotMembArea = 6e-6; %cm^2
        aPR.Units.AreaScaleToM = 1e-4;
        aPR.Units.VoltsScaleToV=1e-3;
        aPR.Units.CondPerUnitAreaToSpercmsqr=1e-3;
        aPR.Units.InjCurrentToAmpsPercmsqr=1e-6;
        aPR.Units.Notes={'all conductances including NMDA and AMPA are in mS/cm^2','All voltages are in mV'};
        aPR.RDSinOhms=1/(aPR.gc*aPR.Units.CondPerUnitAreaToSpercmsqr*aPR.TotMembArea);
        aPR.RDSinMOhms=aPR.RDSinOhms/1e+6;
  %      aPRSS = GetSinglePRSSTonic(1,1000,aPR.Idinj,aPR.Isinj);
  %      aPR.QuasiSS=aPRSS;
        aPR.Units.AreaScaleToM = 1e-4;
        aPR.Units.VoltsScaleToV=1e-3;
        aPR.Units.CondPerUnitAreaToSpercmsqr=1e-3;
        aPR.Units.InjCurrentToAmpsPercmsqr=1e-6;
        aPR.Units.Notes={'all conductances including NMDA and AMPA are in mS/cm^2','All voltages are in mV'};

if nargin > 8
    display('too many variable arguments')
    exit
end
if nargin>=3
    aPR.ENa=varargin{1};
end
if nargin>=4
    aPR.p=varargin{2};
end
if nargin>=5
    aPR.gc=varargin{3};
end
if nargin>=6
    aPR.gKAHP=varargin{4};
end
if nargin>=7
    aPR.gKC=varargin{5};
if nargin>=8
    aPR.Idinj=varargin{6};
end

end