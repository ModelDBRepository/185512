
  function SingPRPlusHWithSyn = SingIntegODE23PRPlusHcurrWithSyn_db(aPRwH,uAmpsPermsecCm2,VsPreSyn,Syn,delay,Tend,VdsOut,VsThresh,SomaInj,termaftevent)

 % 5/24/2015 
 %% What it does:
 % Intgegrates the polarized PR neuron and applies EITHER a current ramp 
 % AMPA or NMDA synaptic currents. This routine uses ODE 23 although ODE45 has also been used. ODE 23 was
 % found to be faster than ODE45 and the solutions were identical to within
 % machine precisions CHECK
 
%%  Dependencies
% GateEquil_db
% GateTimeConst
% alpha and beta functions for h,n,s,c,q
% See PR folder for functions

%% Inputs
% aPR-- is a structure class defined by using InitPR
% uAmpsPermsecCm2-- is the slope of the injected current ramp note units
% delay-- ms of delay before ramp starts
% Tend-- maximum time ode runs
% VdsOut-- the polarization expressed in potential difference outised cell
% from dendrite to soma
% VsThresh--the soma potential at which spike is said to occure and TTFS
% is measured
% SomaInj-- If true injected current ramp goes into soma if False it is
% termafterevent-- to terminate after Vs passes through VsThresh
% applied to dendrite

%%User Inputs for integration parameters. 
%Note an event function is specified
UseDefault=true;
if UseDefault==true
    options = odeset('RelTol',1e-4,'Stats','off','MaxStep',1e-1,'Refine',4,'Events',@events,...
                 'OutputSel',1);
 else
    aRefine=1;
    aRelTol=1e-6;
    aMaxStep=1e-2;
    options = odeset('RelTol',aRelTol,'Stats','on','MaxStep',aMaxStep,'Refine',aRefine,'Events',@events,...
                 'OutputSel',1);   
 end
tic
past=0; %initialize the falg for wether or not a spike has occured

%% Initial state of ODE set prior to this subroutine and stored in aPR.SS
multPRInitY=zeros(1,11);

multPRInitY(1,1)=aPRwH.SS.NumSS(1);
multPRInitY(1,2)=aPRwH.SS.NumSS(2);
multPRInitY(1,3)=aPRwH.SS.NumSS(3);
multPRInitY(1,4)=aPRwH.SS.NumSS(4);
multPRInitY(1,5)=aPRwH.SS.NumSS(5);
multPRInitY(1,6)=aPRwH.SS.NumSS(6);
multPRInitY(1,7)=aPRwH.SS.NumSS(7);
multPRInitY(1,8)=aPRwH.SS.NumSS(8);
multPRInitY(1,9)=0;
multPRInitY(1,10)=0;
multPRInitY(1,11)=aPRwH.SS.NumSS(9);  % CHECK !!!!
%% Reshape to use in MATLAB integrator
multPRInitYCol=zeros(11,1);
multPRInitYCol=reshape(multPRInitY,11,1);

%% Call Ode34 to integrate function defined in PR94NoSyn

[T,YMultcol,te,ye,ie] = ode23(@(t,Y) PR94Syn(t,Y),[0 Tend],multPRInitYCol,options);

% SingPRWithSyn is the structure containing the results of the
% integration with metadat like the filename date and run time.
SingPRPlusHWithSyn.te = te;
SingPRPlusHWithSyn.ye = ye;
SingPRPlusHWithSyn.ie = ie;

SingPRPlusHWithSyn.YMultcol=YMultcol;
SingPRPlusHWithSyn.T = T;


SingPRPlusHWithSyn.etime=toc;
SingPRPlusHWithSyn.datetime=datestr(now);
SingPRPlusHWithSyn.file=mfilename;


SingPRPlusHWithSyn.PR=aPRwH;
SingPRPlusHWithSyn.NumN=1;
SingPRPlusHWithSyn.uAmpsPerCm2=uAmpsPermsecCm2;
SingPRPlusHWithSyn.VdsOut=VdsOut;
SingPRPlusHWithSyn.Tend=Tend;
SingPRPlusHWithSyn.delay=delay;
if isempty(te) 
    SingPRPlusHWithSyn.idxteVs=size(T,1);
else
    SingPRPlusHWithSyn.idxteVs=find(T>=te(1,1),1,'First');
end

%toc
SingPRPlusHWithSyn.etime=toc;
SingPRPlusHWithSyn.datetime=datestr(now);
SingPRPlusHWithSyn.file=mfilename;


function dY = PR94Syn(t,Y)

   Nn=sqrt(size(Y,1));
   Y = reshape(Y,1,11);
if t > delay
if SomaInj==true
   if(Y(1,1) >= VsThresh | past == 1)
        Isinj =aPRwH.Isinj;
        past = 1;
   else
        Isinj=aPRwH.Isinj+heaviside(t-delay)*uAmpsPermsecCm2*(t-delay);
   end
   Idinj=aPRwH.Idinj;
else
    if Y(1,1) >= VsThresh
       Idinj =aPRwH.Idinj;
    else
        Idinj=aPRwH.Idinj+heaviside(t-delay)*uAmpsPermsecCm2*(t-delay);
    end
   Isinj=aPRwH.Isinj;
end
else %t less than delay
    Isinj=aPRwH.Isinj;
    Idinj=aPRwH.Idinj;
end
    Cm=aPRwH.Cm;
    gL=aPRwH.gL;
    gNa=aPRwH.gNa;
    gKDR=aPRwH.gKDR;
    gKC=aPRwH.gKC;
    gKAHP=aPRwH.gKAHP;
    gCa=aPRwH.gCa;
    ENa=aPRwH.ENa;            %CHECK not sure this is good idea to double the # variable inside function to be integrated
    Ek=aPRwH.Ek;
    EL=aPRwH.EL;
    ECa=aPRwH.ECa;
    p=aPRwH.p;
    gc=aPRwH.gc;
    WRT=aPRwH.WRT;
    Vsyn=aPRwH.Vsyn;
    MaxS=aPRwH.MaxS;
    gNMDA=Syn.gNMDA;
    gAMPA=Syn.gAMPA;
    HcurrWRT=aPRwH.HcurrLippertWRT;
    h_Vhalf=aPRwH.h_Vhalf;
    VsPreAmp=VsPreSyn.Amp;
    t_initexc=VsPreSyn.t_initexc;
    t_spkdur=VsPreSyn.t_spkdur;
if Y(9) > MaxS
      Y(9)=MaxS;
end;
expNMDA=1./(1+0.28*exp(-0.062*(Y(2)-60)));
INMDA=gNMDA*Y(2).*Y(9).*expNMDA-gNMDA*Vsyn*Y(9).*expNMDA;
IAMPA=gAMPA*Y(10).*Y(2)-gAMPA*Vsyn*Y(10);
VsPre=VsPreAmp*heaviside(t-t_initexc)*heaviside(t_initexc+t_spkdur-t);
dY = zeros(1,10);    % a column vector

dY(1) = (1/Cm)*(-gL*(Y(1)-EL)-gNa*MInfPR94(Y(1),WRT).*Y(4).*(Y(1)-ENa)-gKDR*Y(5).*(Y(1)-Ek)...
        +(gc/p)*(Y(2)-Y(1)+VdsOut)+Isinj/p);
dY(2) = (1/Cm)*(-gL*(Y(2)-EL)-gCa*(Y(2)-ECa).*Y(6).^2-gKAHP*Y(8).*(Y(2)-Ek)-gKC*Y(7).*Chi(Y(3)).*(Y(2)-Ek)...
            -  aPRwH.gh*Y(11)*(Y(2)-aPRwH.Eh)+gc/(1-p)*(Y(1)-Y(2)-VdsOut)+Idinj/(1-p)-(INMDA+IAMPA)/(1-p));

dY(3) = -0.13*gCa*(Y(2)-ECa).*Y(6).^2-0.075*Y(3);
dY(4) = (GateEquil_db(alphah_db(Y(1),WRT),betah_db(Y(1),WRT))-Y(4))./GateTimeCnst_db(alphah_db(Y(1),WRT),betah_db(Y(1),WRT));
dY(5) = (GateEquil_db(alphan_db(Y(1),WRT),betan_db(Y(1),WRT))-Y(5))./GateTimeCnst_db(alphan_db(Y(1),WRT),betan_db(Y(1),WRT));
dY(6) = (GateEquil_db(alphas_db(Y(2),WRT),betas_db(Y(2),WRT))-Y(6))./GateTimeCnst_db(alphas_db(Y(2),WRT),betas_db(Y(2),WRT));
dY(7) = (GateEquil_db(alphac_db (Y(2),WRT),betac_db(Y(2),WRT))-Y(7))./GateTimeCnst_db(alphac_db(Y(2),WRT),betac_db(Y(2),WRT));
dY(8) = (GateEquil_db(alphaq_db(Y(3)),betaq_db)-Y(8))./GateTimeCnst_db(alphaq_db(Y(3)),betaq_db);
dY(9) = SynSumNMDA(VsPre) - Y(9)/150;
dY(10) = SynSumAMPA(VsPre) - Y(10)/2; 
dY(11) = (GateEquil_db(alpha_idB(Y(2),h_Vhalf,HcurrWRT),beta_idB(Y(2),h_Vhalf,HcurrWRT))-Y(11))./GateTimeCnst_db(alpha_idB(Y(2),h_Vhalf,HcurrWRT),beta_idB(Y(2),h_Vhalf,HcurrWRT));
dY = reshape(dY,11,1);
end
function SumNMDA = SynSumNMDA(VsPre)
    SumNMDA = 0;
    if VsPre > 10 
        SumNMDA=SumNMDA+1;
    end
end

function SumAMPA = SynSumAMPA(VsPre)
    SumAMPA = 0;
    if VsPre > 20
        SumAMPA=SumAMPA+1;
    end
end
function CaSatChi= Chi(Ca)
    CaSatChi = min(Ca/250,1);
end

function MInfsqr = MInfPR94(Vs,WRT) 
    alp=alpham_db(Vs,WRT);
    bet=betam_db(Vs,WRT);
    MInfsqr = GateEquil_db(alp,bet).^2;
end

function [value,isterminal,direction] = events(t,y)
% Locate the time when potential passes through zero in a 
% decreasing direction and stop integration.
value = y(1)-VsThresh;     % Detect Soma = 10
if termaftevent==true
    isterminal = 1;   % Stop the integration
else
    isterminal = 0;   % Keep going
end
direction = 1;   % positive direction only
end
end
