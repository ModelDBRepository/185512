%% Figure 3. 
% Weak and Intermediate polarization 
% With Linear Fit. 
% Assumes four panels each with a different _pair_ of M and E_{K}
% Default: Ek={-25,-45} mV X M = {0.3,0.8} \muA/(cm^{2} s)
% Ramp current injection protocol
% RIR July 5, 2015

%% Dependencies
% # IniPR_db
% # NumerEquilPR_db
% # SingIntegODE23PRWithInjCurr_db

%% User Supplied Inputs
% # Array of V_{ds}^{out} to plot
%   MaxVdsOut, DeltaVdsOut, NumVdsOut 
% 


 
MaxVdsOut=10; %Default: 10 mV
DeltaVdsOut=0.5; %Default: 0.5 mV
NumVdsOut=3;   %Default: 51
Ms=0.8;
Eks=-25;
guessVsVd=[0,0];
SomaInj=true;
%VsThresh=30;
VsThresh=10;
%delay =500; 
delay=40; %prior to this model neuron should be at rest with VdsOut. Provides a good sanity check clear before and after, just remember to subtract delay from spike time to get TTFS
Tend=17500;
termaftevent=false;
Isinj=-0.5;
smp=[15,47];
f1=figure();
for i=1:size(Ms,2)
    tmpM=Ms(i);
    tmpuAmpsPermsecCm2=tmpM*1e-3;
    for j=1:size(Eks,2)
       tmpEk=Eks(j); 
       tmpPR=IniPR_db(Isinj,tmpEk);
       for k=1:NumVdsOut
           tmpVdsOut=MaxVdsOut-(k-1)*DeltaVdsOut;
           [tmpnumSSPR,diffProjFullEq,Jacob,eigJacob,nzeig] = NumerEquilPR_db(tmpPR,guessVsVd,tmpVdsOut); 
           tmpPR.SS.NumSS=tmpnumSSPR;
           tmpSingPR = SingIntegODE23PRWithInjCurr_db(tmpPR,tmpuAmpsPermsecCm2,delay,Tend,tmpVdsOut,VsThresh,SomaInj,termaftevent);
           subsample=100;
           extr=1:subsample:size(tmpSingPR.T,1);
           VsVdTrace(i,j).T(k).T=tmpSingPR.T(extr);
           VsVdTrace(i,j).Vs(k).Vs=tmpSingPR.YMultcol(extr,1);
           VsVdTrace(i,j).Vd(k).Vd=tmpSingPR.YMultcol(extr,2);
           VsVdTrace(i,j).VdsOut(k)=tmpVdsOut;
           VsVdTrace(i,j).TTFS(k)=tmpSingPR.te(1)-tmpSingPR.delay;
           VsVdTrace(i,j).NumSpks(k)=size(tmpSingPR.te,1);
           VsVdTrace(i,j).FirstToLastSpike(k)=tmpSingPR.te(size(tmpSingPR.te,1))-tmpSingPR.te(1);
           VsVdTrace(i,j).M=tmpM;
           VsVdTrace(i,j).Ek=tmpEk;
           k
       end
%        subplot(2,1,1)
%        def=linspace(1,NumVdsOut,NumVdsOut);
%  %      plot(VsVdTrace(i,j).VdsOut(nosubidx),VsVdTrace(i,j).TTFS(nosubidx),'ok','MarkerFaceColor','k','MarkerSize',2)
%        plot(VsVdTrace(i,j).VdsOut(extr(1,:)),VsVdTrace(i,j).TTFS(extr(1,:)),'ok','MarkerFaceColor','k','MarkerSize',2)
%        title({['m = ', num2str(tmpM),' \muA/(cm^{2} s)',' Ek= ',num2str(tmpEk),' mV',' VsThresh = ',num2str(VsThresh),' mV']},'FontSize',10)
%        subplot(2,1,2)
%        if mod(k,10)==0
%             plot(VsVdTrace(i,j).T(smp(1)).T,VsVdTrace(i,j).Vs(smp(1)).Vs,'ok','MarkerFaceColor','w','MarkerSize',2,'MarkerEdgeColor','w')
%             hold on;
%        end
%        ax3=gca();
%        set(ax3,'FontSize',8)
    end
end

f2=figure();

subplot(2,1,1)
plot(VsVdTrace(1,1).VdsOut(1,:),VsVdTrace(1,1).NumSpks(1,:),'sk')
title({['m = ', num2str(tmpM),' \muA/(cm^{2} s)',' Ek= ',num2str(tmpEk),' mV',' VsThresh = ',num2str(VsThresh),' mV']},'FontSize',10)
subplot(2,1,2)
plot(VsVdTrace(1,1).VdsOut(1,:),VsVdTrace(1,1).FirstToLastSpike(1,:),'-dk')

f3=figure();
title({['m = ', num2str(tmpM),' \muA/(cm^{2} s)',' Ek= ',num2str(tmpEk),' mV',' VsThresh = ',num2str(VsThresh),' mV']},'FontSize',10)
for k=1:size(VsVdTrace(1,1).VdsOut,2)
    if VsVdTrace(1,1).NumSpks(1,k)==0
        sym='ok';
        col='k';
    elseif VsVdTrace(1,1).NumSpks(1,k)==1
        sym='sk';
        col='None';
    elseif VsVdTrace(1,1).NumSpks(1,k)==2
        sym='^k';
        col='k';
    elseif VsVdTrace(1,1).NumSpks(1,k)==3
        sym='dk';
        col='None';
    elseif VsVdTrace(1,1).NumSpks(1,k) < 10 && VsVdTrace(1,1).NumSpks(1,k) > 3
        sym='+k';
        col='k';
    else
        sym='*k';
        col='None';     
    end
    plot(VsVdTrace(1,1).VdsOut(1,k),VsVdTrace(1,1).TTFS(1,k),sym,'MarkerFaceColor',col,'MarkerSize',6)
    hold on;
end
