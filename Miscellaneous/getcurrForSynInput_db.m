function [curr pot Spikes] = getcurrForSynInput_db(aSynPR,Syn,VdsOut,uAmpsPermsecCm2,delay)

%% RIR June 7, 2015
%% What does it do?



SomaInj=true;

NumT=size(aSynPR.YMultcol,1);
curr.NumT=NumT;
curr.IsLeak=zeros(NumT,1);
curr.IsNa=zeros(NumT,1);
curr.IsKdr=zeros(NumT,1);
curr.IsVdVs=zeros(NumT,1);
curr.Isinj=zeros(NumT,1);
curr.IdLeak=zeros(NumT,1);
curr.IdCa=zeros(NumT,1);
curr.IdKAHP=zeros(NumT,1);
curr.IdKC=zeros(NumT,1);
curr.IdVsVd=zeros(NumT,1);
curr.Idinj=zeros(NumT,1);
curr.T=zeros(NumT,1);
curr.IAMPA=zeros(NumT,1);
curr.INMDA=zeros(NumT,1);
curr.TotalSoma=zeros(NumT,1);
curr.TotalDend=zeros(NumT,1);
pot.Vs=zeros(NumT,1);
pot.Vd=zeros(NumT,1);
pot.Ca = zeros(NumT,1);

curr.T=aSynPR.T;
curr.VdsOut=VdsOut;
%curr.SomaInj=SomaInj;
curr.IsinjBase=aSynPR.PR.Isinj;
curr.IdinjBase=aSynPR.PR.Idinj;
curr.uAmpsPermsecCm2=uAmpsPermsecCm2;
curr.delay=delay;
curr.Ek=aSynPR.PR.Ek;

pot.Vs=aSynPR.YMultcol(:,1);
pot.Vd=aSynPR.YMultcol(:,2);
pot.Ca=aSynPR.YMultcol(:,3);
pot.h=aSynPR.YMultcol(:,4);
pot.n=aSynPR.YMultcol(:,5);
pot.s=aSynPR.YMultcol(:,6);
pot.c=aSynPR.YMultcol(:,7);
pot.q=aSynPR.YMultcol(:,8);
pot.S=aSynPR.YMultcol(:,9);
pot.W=aSynPR.YMultcol(:,10);

curr.IsLeak=-aSynPR.PR.gL*(aSynPR.YMultcol(:,1)-aSynPR.PR.EL);
curr.IsNa=-aSynPR.PR.gNa*MInfPR94(aSynPR.YMultcol(:,1),aSynPR.PR.WRT).*aSynPR.YMultcol(:,4).*(aSynPR.YMultcol(:,1)-aSynPR.PR.ENa);
curr.IsKdr=-aSynPR.PR.gKDR*aSynPR.YMultcol(:,5).*(aSynPR.YMultcol(:,1)-aSynPR.PR.Ek);
%Belwo not correct for ephaptic
curr.IsVdVs=(aSynPR.PR.gc/aSynPR.PR.p)*(aSynPR.YMultcol(:,2)-aSynPR.YMultcol(:,1)+VdsOut);
expNMDA=1./(1+0.28*exp(-0.062*(aSynPR.YMultcol(:,2)-60)));
curr.INMDA=Syn.gNMDA*aSynPR.YMultcol(:,2).*aSynPR.YMultcol(:,9).*expNMDA-Syn.gNMDA*aSynPR.PR.Vsyn*aSynPR.YMultcol(:,9).*expNMDA;
curr.IAMPA=Syn.gAMPA*aSynPR.YMultcol(:,10).*aSynPR.YMultcol(:,2)-Syn.gAMPA*aSynPR.PR.Vsyn*aSynPR.YMultcol(:,10);

curr.IdLeak=-aSynPR.PR.gL*(aSynPR.YMultcol(:,2)-aSynPR.PR.EL);
curr.IdCa=-aSynPR.PR.gCa*(aSynPR.YMultcol(:,2)-aSynPR.PR.ECa).*aSynPR.YMultcol(:,6).^2;
curr.IdKAHP=-aSynPR.PR.gKAHP*aSynPR.YMultcol(:,8).*(aSynPR.YMultcol(:,2)-aSynPR.PR.Ek);
curr.IdKC=-aSynPR.PR.gKC*aSynPR.YMultcol(:,7).*Chi(aSynPR.YMultcol(:,3)).*(aSynPR.YMultcol(:,2)-aSynPR.PR.Ek);
%Belwo not correct for ephaptic
curr.IdVsVd=aSynPR.PR.gc/(1-aSynPR.PR.p)*(aSynPR.YMultcol(:,1)-aSynPR.YMultcol(:,2)-VdsOut);
curr.q=aSynPR.YMultcol(:,8);

pastdelay=false;
injidx=0;
if isempty(aSynPR.te)
       Spikes=false;
else
    Spikes=true;
end
for i=1:NumT
t=aSynPR.T(i,1);
if SomaInj==true && Spikes
    if t < aSynPR.te(1,1) && t > delay
        if pastdelay==true
            
        else
            injidx=i;
            pastdelay=true;
        end
        curr.Isinj(i)=(aSynPR.PR.Isinj+uAmpsPermsecCm2*(t-delay))/aSynPR.PR.p;
    else
        curr.Isinj(i)=aSynPR.PR.Isinj/aSynPR.PR.p;
    end
   %curr.Idinj=aSynPR.PR.Idinj.*ones(size(curr.Isinj,1),1);
   curr.Idinj(i)=aSynPR.PR.Idinj;
elseif SomaInj==false && Spikes==true;
    if t < aSynPR.te(1,1) && t > delay
        curr.Idinj(i)=(aSynPR.PR.Idinj+uAmpsPermsecCm2*(t-delay))/(1-aSynPR.PR.p);
    else
        injidx=i;
        curr.Idinj(i)=aSynPR.PR.Idinj/(1-aSynPR.PR.p);
    end
   %curr.Isinj=aSynPR.PR.Isinj.*ones(size(curr.Idinj,1),1);
   curr.Isinj(i)=aSynPR.PR.Isinj;
end
end
aSynPR.idxteVs=size(curr.T,1);


%%%% For tomorrow need to change denominator.

injidx=injidx+1;
aSynPR.idxteVs=aSynPR.idxteVs-1;
%strtidx=injidx;
strtidx=1;
curr.JustSomaInj=curr.Isinj(strtidx:aSynPR.idxteVs)-aSynPR.PR.Isinj;
curr.TotalSoma = curr.IsLeak(strtidx:aSynPR.idxteVs)+curr.IsNa(strtidx:aSynPR.idxteVs)+curr.IsKdr(strtidx:aSynPR.idxteVs)+curr.IsVdVs(strtidx:aSynPR.idxteVs)+curr.Isinj(strtidx:aSynPR.idxteVs);
curr.TotalDend = curr.IdLeak(strtidx:aSynPR.idxteVs)+curr.IdCa(strtidx:aSynPR.idxteVs)+curr.IdKAHP(strtidx:aSynPR.idxteVs)+curr.IdKC(strtidx:aSynPR.idxteVs)+curr.IdVsVd(strtidx:aSynPR.idxteVs)+curr.Idinj(strtidx:aSynPR.idxteVs);
curr.TotalActiveDend=curr.IdCa(strtidx:aSynPR.idxteVs)+curr.IdKAHP(strtidx:aSynPR.idxteVs)+curr.IdKC(strtidx:aSynPR.idxteVs);
curr.TotalActiveSoma=curr.IsNa(strtidx:aSynPR.idxteVs)+curr.IsKdr(strtidx:aSynPR.idxteVs);
curr.TotalSyn=curr.INMDA(strtidx:aSynPR.idxteVs)+curr.IAMPA(strtidx:aSynPR.idxteVs);
curr.IsVdVsInjIs=curr.IsVdVs(strtidx:aSynPR.idxteVs)./curr.JustSomaInj;
curr.IdVsVdInjIs=curr.IdVsVd(strtidx:aSynPR.idxteVs)./curr.JustSomaInj;
curr.begidx=strtidx;
curr.endidx=aSynPR.idxteVs;

curr.TotalInOut = curr.TotalSoma+curr.TotalDend;
curr.FracKAHPDend=curr.IdKAHP(strtidx:aSynPR.idxteVs)./curr.TotalDend;
curr.FracKAHPSoma=curr.IdKAHP(strtidx:aSynPR.idxteVs)./curr.TotalSoma;
curr.FracKAHPInjSoma=curr.IdKAHP(strtidx:aSynPR.idxteVs)./curr.JustSomaInj;
curr.FracKAHPJustInjSoma=curr.IdKAHP(strtidx:aSynPR.idxteVs)./curr.JustSomaInj;
curr.FracKAHPVsVd=curr.IdKAHP(strtidx:aSynPR.idxteVs)./curr.JustSomaInj;



function CaSatChi= Chi(Ca)
    CaSatChi = min(Ca/250,1);
end

function MInfsqr = MInfPR94(Vs,WRT) 
    alp=alpham_db(Vs,WRT);
    bet=betam_db(Vs,WRT);
    MInfsqr = GateEquil_db(alp,bet).^2;
end
end

