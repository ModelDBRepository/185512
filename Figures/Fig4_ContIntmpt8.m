%% Figure # in paper  
% RIR July 15, 2015

%% What it does
% code used to create figs 4 and 5 in paper
% Using injected ramp protocol for pairs of M and Ek it plots TTFS versus
% VdsOut as in Fig. 3 . However, once Vs reaches 30 mV the current stays at
% a constant level and we keep integrating. We do this to get a feel for
% the different polarization dependent spiking dynamics.


tic
%% Run Paramters
% Define a spike by Vs reaching a certain threshold VsThresh

VsThresh=30;
%VsThresh=10; %The TTFS should not be very sensitive to VsThresh. Except in
%the case of very strongly negative polarizations or if the neuron goes
%into a bursting like mode
delay=40; % This is extra extra precaution just inn case its not quite at its equilibrium. This has never really proved neceessary
            % Remember to get TTFS we need to subtract dealy from event
            % time of Vs passing through VsThresh
%Tend=7500; %These are the maximum time the rampinjection is applied the
% as long as IsTerm=true then the computation at that polarization will
% stop when Vs passes through VsThresh
%Tend=12000;
Tend=10000;
%Tend=250; %m/s much quicker than ramp injected currents at least for the M we had been using. Actually most are under 100 ms
SomaInj=true; %Optionally the current ramp may be injected in the dendrite for gc=2.1 and p=0.5 there is not much qualitative difference
                % This does not matter since for this syanptic stimulus we
                % will set M=0;
termafterevent=true; % Stops after Vs first rises up and passes through VsThresh
DoNotTermafterEvent=false;
%% Initialize PR neuron

Idinj=0; % Nothing extra is injected in the dendrite besides the synaptic current
%Isinj=-0.3; %uA/cm^2This should be selected to match known spike intervals times with injected current. -0.5 is the one used in PR and Grahm, Cutridis and Hunter in Associative model 
Isinj=-0.5;
% chapter in Hippocampal Microcircuits use -0.5 as. The origins of -0.3 is somehat of a mistake I used that when I also had FSI going. Would like to run it again sometime with Is=-0.5
% fits several other computational models using synaptic input and the PR model (Need References) 
M=0.8;
tmpuAmpsPermsecCm2=M*1e-3; %per msec
%NumtmpVdsOut=33;
%NumtmpVdsOut=15; %from -12.25 this goes to -16
%NumtmpVdsOut=75; %from -12.25 this goes to -16
NumVdsOut=15;
VdsTTFS1=zeros(NumVdsOut,2);
%StarttmpVdsOut=-4;
%StarttmpVdsOut=-12.25;
StarttmpVdsOut=0;
%StarttmpVdsOut=-16.25;
SteptmpVdsOut=-1;
%SteptmpVdsOut=-0.075;
% MinVdsOyt=strttmpVdsOut+NumtmpVdsOut*SteptmpVdsOut 
MintmpVdsOut=StarttmpVdsOut+NumVdsOut*SteptmpVdsOut;
% error VdsTTFSEkgAMPA=zeros(NumtmpVdsOut,2); found May 29, 2015
%% Integration Parameters
RelTol=1e-4;
AbsTol=1e-6;
MaxStep=1e-1;
IntegParams=[RelTol,AbsTol,MaxStep];

guessVsVd(1)=-10; % To initiate finding the equilibrium. -10, -10 is as good a place as any not real sensitive
guessVsVd(2)=-10;



%% Define 4 sets of $E_{K}$ and Initialize 4 PR neurons values

Ek1=-25.0; %mV
Ek2=-45.0; %mV


aPR1=IniPR_db(Isinj,Ek1);
aPR2=IniPR_db(Isinj,Ek2);


%% collection of symbols
symbag={'-sk','-ok'};

%% Loop over range of polarizations.

% User specified 3 polarizations to plot
samplevds=[-2,-7,-15];
tafter=10000; %ms Also means that for isoloated spikes if it is actually periodic it has a frequency less than 0.1 Hs
MaxPulseLength=20; %ms this means that the fastest frequency of bursts is 1/MaxPulseLength for 20 ms it is 50 Hz
for i = 1:NumVdsOut
i
tmpVdsOut = StarttmpVdsOut+(SteptmpVdsOut*i);

[numSSPR1,diffProjFullEq,Jacob1,eigJacob1,nzeig] = NumerEquilPR_db(aPR1,guessVsVd,tmpVdsOut);

aPR1.SS.NumSS=numSSPR1;
aPR1.SS.eig = eig(Jacob1);
if aPR1.SS.eig < 0 %check if all eigenvalues less than zero
    stable1=1;
    if tmpVdsOut==samplevds(1) || tmpVdsOut==samplevds(2) || tmpVdsOut==samplevds(3)
          [SingPRWithInjCurrbefore, SingPRWithInjCurrafter] =  SingIntegODE23PRWithInjCurrWintegOptAftCurrCnst_db(aPR1,tmpuAmpsPermsecCm2,delay,Tend,tmpVdsOut,VsThresh,SomaInj,DoNotTermafterEvent,IntegParams,tafter);
    else
          [SingPRWithInjCurrbefore, SingPRWithInjCurrafter] =  SingIntegODE23PRWithInjCurrWintegOptAftCurrCnst_db(aPR1,tmpuAmpsPermsecCm2,delay,Tend,tmpVdsOut,VsThresh,SomaInj,DoNotTermafterEvent,IntegParams,tafter);
    end
   subsample=10;
   SingPRWithInjCurr1.T=[SingPRWithInjCurrbefore.T',SingPRWithInjCurrafter.T']';
   SingPRWithInjCurr1.YMultcol=[SingPRWithInjCurrbefore.YMultcol',SingPRWithInjCurrafter.YMultcol']';
   extr=1:subsample:size(SingPRWithInjCurr1.T,1);
   VsVdTrace1.T(i).T=SingPRWithInjCurr1.T(extr);
   VsVdTrace1.Vs(i).Vs=SingPRWithInjCurr1.YMultcol(extr,1);
   VsVdTrace1.Vd(i).Vd=SingPRWithInjCurr1.YMultcol(extr,2);
   VsVdTrace1.VdsOut(i)=tmpVdsOut;
   VsVdTrace1.TTFS(i)=SingPRWithInjCurrbefore.te(1)-SingPRWithInjCurrbefore.delay;
   VsVdTrace1.NumSpks(i)=size(SingPRWithInjCurrbefore.te,1)+size(SingPRWithInjCurrafter.te,1);  %DOes passing throug Vs on the way up before  mean that it will not be conte in after? Assume yes
   VsVsTrace1.DiffSpk(i).diff=diff([SingPRWithInjCurrbefore.te',SingPRWithInjCurrafter.te'])';
   VsVdTrace1.DiffSpkGTMaxPulse(i)=size(find(diff([SingPRWithInjCurrbefore.te',SingPRWithInjCurrafter.te'])'>MaxPulseLength),1);
   VsVdTrace1.NumPulses(i)=size(find(diff([SingPRWithInjCurrbefore.te',SingPRWithInjCurrafter.te'])'>MaxPulseLength),1)+1;
   VsVdTrace1.DiffSpkLTMaxPulse(i)=size(find(diff([SingPRWithInjCurrbefore.te',SingPRWithInjCurrafter.te'])'<=MaxPulseLength),1);
   VsVdTrace1.NumSpiklets(i)=size(find(diff([SingPRWithInjCurrbefore.te',SingPRWithInjCurrafter.te'])'<=MaxPulseLength),1);
   VsVdTrace1.NumSpikletsPerPulse(i)=VsVdTrace1.NumSpiklets/VsVdTrace1.NumPulses;
   if isempty(SingPRWithInjCurrafter.te)
       VsVdTrace1.FirstToLastSpike(i)=NaN;
   else
       VsVdTrace1.FirstToLastSpike(i)=SingPRWithInjCurrafter.te(size(SingPRWithInjCurrafter.te,1))-SingPRWithInjCurrbefore.te(1);
   end
   
  % VsVdTrace1.FirstToLastSpike(i)=SingPRWithInjCurrafter.te(size(SingPRWithInjCurrafter.te,1))-SingPRWithInjCurrbefore.te(1);
   VsVdTrace1.M=M;
   VsVdTrace1.Ek=Ek1;
   % SingPRWithSyn.datetime
else
    stable1=0;
end
SomaInj=true;
if stable1 == 1
    if isempty(SingPRWithInjCurrbefore.te(1)) 
        VdsTTFS1(i,1)=tmpVdsOut;
        VdsTTFS1(i,2)=NaN; %NaN didn't spike
    else
        VdsTTFS1(i,1)=tmpVdsOut;
        VdsTTFS1(i,2)=SingPRWithInjCurrbefore.te(1)-SingPRWithInjCurrbefore.delay; 
    end
else
    VdsTTFS1(i,1)=tmpVdsOut; 
    VdsTTFS1(i,2)=-1; %unstable
end
%% 2nd Neuron
[numSSPR2,diffProjFullEq,Jacob2,eigJacob2,nzeig] = NumerEquilPR_db(aPR2,guessVsVd,tmpVdsOut);

aPR2.SS.NumSS=numSSPR2;
aPR2.SS.eig = eig(Jacob2);
if aPR2.SS.eig < 0 %check if all eigenvalues less than zero
    stable2=1;
    if tmpVdsOut==samplevds(1) || tmpVdsOut==samplevds(2) || tmpVdsOut==samplevds(3)
          [SingPRWithInjCurrbefore2, SingPRWithInjCurrafter2] =  SingIntegODE23PRWithInjCurrWintegOptAftCurrCnst_db(aPR2,tmpuAmpsPermsecCm2,delay,Tend,tmpVdsOut,VsThresh,SomaInj,DoNotTermafterEvent,IntegParams,tafter);
    else
          [SingPRWithInjCurrbefore2, SingPRWithInjCurrafter2] =  SingIntegODE23PRWithInjCurrWintegOptAftCurrCnst_db(aPR2,tmpuAmpsPermsecCm2,delay,Tend,tmpVdsOut,VsThresh,SomaInj,DoNotTermafterEvent,IntegParams,tafter);
    end
   subsample=10;
   SingPRWithInjCurr2.T=[SingPRWithInjCurrbefore2.T',SingPRWithInjCurrafter2.T']';
   SingPRWithInjCurr2.YMultcol=[SingPRWithInjCurrbefore2.YMultcol',SingPRWithInjCurrafter2.YMultcol']';
   extr=1:subsample:size(SingPRWithInjCurr2.T,1);
   VsVdTrace2.T(i).T=SingPRWithInjCurr2.T(extr(1,:));
   VsVdTrace2.Vs(i).Vs=SingPRWithInjCurr2.YMultcol(extr(1,:),1);
   VsVdTrace2.Vd(i).Vd=SingPRWithInjCurr2.YMultcol(extr(1,:),2);
   VsVdTrace2.VdsOut(i)=tmpVdsOut;
   VsVdTrace2.TTFS(i)=SingPRWithInjCurrbefore2.te(1)-SingPRWithInjCurrbefore2.delay;
   VsVdTrace2.NumSpks(i)=size(SingPRWithInjCurrbefore2.te,1)+size(SingPRWithInjCurrafter2.te,1);  %DOes passing throug Vs on the way up before  mean that it will not be conte in after? Assume yes
   VsVdTrace2.DiffSpk(i).diff=diff([SingPRWithInjCurrbefore2.te',SingPRWithInjCurrafter2.te'])';
   VsVdTrace2.DiffSpkGTMaxPulse(i)=size(find(diff([SingPRWithInjCurrbefore2.te',SingPRWithInjCurrafter2.te'])'>MaxPulseLength),1);
   VsVdTrace2.NumPulses(i)=size(find(diff([SingPRWithInjCurrbefore2.te',SingPRWithInjCurrafter2.te'])'>MaxPulseLength),1)+1;
   VsVdTrace2.DiffSpkLTMaxPulse(i)=size(find(diff([SingPRWithInjCurrbefore2.te',SingPRWithInjCurrafter2.te'])'<=MaxPulseLength),1);
   VsVdTrace2.NumSpiklets(i)=size(find(diff([SingPRWithInjCurrbefore2.te',SingPRWithInjCurrafter2.te'])'<=MaxPulseLength),1);
   VsVdTrace2.NumSpikletsPerPulse(i)=VsVdTrace2.NumSpiklets(i)/VsVdTrace2.NumPulses(i);
   
   if isempty(SingPRWithInjCurrafter2.te)
       VsVdTrace2.FirstToLastSpike(i)=NaN;
   else
       VsVdTrace2.FirstToLastSpike(i)=SingPRWithInjCurrafter2.te(size(SingPRWithInjCurrafter2.te,1))-SingPRWithInjCurrbefore2.te(1);
   end
   VsVdTrace2.M=M;
   VsVdTrace2.Ek=Ek2;
    % SingPRWithSyn.datetime
else
    stable2=0;
end
SomaInj=true;
if stable2 == 1
    if isempty(SingPRWithInjCurrbefore2.te) 
        VdsTTFS2(i,1)=tmpVdsOut;
        VdsTTFS2(i,2)=NaN; %NaN didn't spike
    else
        VdsTTFS2(i,1)=tmpVdsOut;
        VdsTTFS2(i,2)=SingPRWithInjCurrbefore2.te(1)-SingPRWithInjCurrbefore2.delay; 
    end
else
    VdsTTFS2(i,1)=tmpVdsOut; 
    VdsTTFS2(i,2)=-1; %unstable
end


end % end loop over range of polarization

f1=figure();
subplot(6,3,[1,4,7,10,13,16])

for k=1:size(VsVdTrace1.VdsOut,2)
    if VsVdTrace1.NumPulses(k)==1
        %single isolated spike or burst
        col='None';
    elseif VsVdTrace1.NumPulses(k)>1
        col='k';
    else
        col='None';
    end
    
    if round(VsVdTrace1.NumSpikletsPerPulse(k))==1
        sym='ok';
    elseif round(VsVdTrace1.NumSpikletsPerPulse(k))==2
        sym='^k';
    elseif round(VsVdTrace1.NumSpikletsPerPulse(k))>=3
        sym='sk';
    elseif  round(VsVdTrace1.NumSpikletsPerPulse(k))==0
        %it there are no spikes which should be impossibe  both since we are using ramp injection protocol and also since there is 1 pulse
        sym='dk';
    else
        sym='*k';
    end
    plot(VsVdTrace1.VdsOut(1,k),VsVdTrace1.TTFS(1,k),sym,'MarkerFaceColor',col,'MarkerSize',6)
    hold on;
end
xlabel('V_{ds}^{out} (mV)')
ylabel('TTFS (ms)')
% for k=1:size(VsVdTrace1.VdsOut,2)
%     if VsVdTrace1.NumSpks(1,k)==0
%         sym='ok';
%         col='k';
%     elseif VsVdTrace1.NumSpks(1,k)==1
%         sym='sk';
%         col='None';
%     elseif VsVdTrace1.NumSpks(1,k)==2
%         sym='^k';
%         col='k';
%     elseif VsVdTrace1.NumSpks(1,k)==3
%         sym='dk';
%         col='None';
%     elseif VsVdTrace1.NumSpks(1,k) < 10 && VsVdTrace1.NumSpks(1,k) > 3
%         sym='+k';
%         col='k';
%     else
%         sym='*k';
%         col='None';     
%     end
%     plot(VsVdTrace1.VdsOut(1,k),VsVdTrace1.TTFS(1,k),sym,'MarkerFaceColor',col,'MarkerSize',6)
%     hold on;
% end

for k=1:size(VsVdTrace2.VdsOut,2)
    if VsVdTrace2.NumPulses(k)==1
        %single isolated spike or burst
        col='None';
    elseif VsVdTrace2.NumPulses(k)>1
        col='k';
    else
        col='None';
    end
    
    if round(VsVdTrace2.NumSpikletsPerPulse(k))==1
        sym='ok';
    elseif round(VsVdTrace2.NumSpikletsPerPulse(k))==2
        sym='^k';
    elseif round(VsVdTrace2.NumSpikletsPerPulse(k))>=3
        sym='sk';
    elseif  round(VsVdTrace2.NumSpikletsPerPulse(k))==0
        %it there are no spikes which should be impossibe  both since we are using ramp injection protocol and also since there is 1 pulse
        sym='dk';
    else
        sym='*k';
    end
    
    plot(VsVdTrace2.VdsOut(1,k),VsVdTrace2.TTFS(1,k),sym,'MarkerFaceColor',col,'MarkerSize',6)
    hold on;
end

% for k=1:size(VsVdTrace2.VdsOut,2)
%     if VsVdTrace2.FirstToLastSpike(1,k) > MinPulseInterval
%             col='k';
%             multfactor=floor(tafter/VsVdTrace2.FirstToLastSpike(1,k))+1;
%     else
%             col='none';
%             multfactor=1;
%     end
%     if VsVdTrace2.NumSpks(1,k)==0
%         sym='+k';       
%         col='None';
%     elseif VsVdTrace2.NumSpks(1,k)==1
%         sym='ok';
%         col='None';
%     elseif VsVdTrace2.NumSpks(1,k)==2
%         sym='^k';
%         
%     elseif VsVdTrace2.NumSpks(1,k)==3
%         sym='sk';
%        
%     elseif VsVdTrace2.NumSpks(1,k) < 10 && VsVdTrace2.NumSpks(1,k) > 3
%         sym='+k';
%         col='k';
%     else
%         sym='*k';
%         col='None';     
%     end
%     plot(VsVdTrace2.VdsOut(1,k),VsVdTrace2.TTFS(1,k),sym,'MarkerFaceColor',col,'MarkerSize',6)
%     hold on;
% end




% subplot(3,2,[1,3,5])
% plot(VdsTTFS1(:,1),VdsTTFS1(:,2),symbag{1},VdsTTFS2(:,1),VdsTTFS2(:,2),symbag{2})
% legend(['E_{K}= ',num2str(Ek1),' (mV) M = ', num2str(M),' mS/cm^{2}'],['E_{K}= ',num2str(Ek2),' (mV) M = ', num2str(M),' mS/cm^{2}'])
% title(['I_{s}= ',num2str(Isinj),' muA/cm^{2}'])
% xlim([-15,4]);
%% Begin computation and plotting of soma potential traces



subplot(6,3,2)
%Ek1 Vds1
tmpVds1idx=find(VsVdTrace1.VdsOut(1,:)==samplevds(1));
tmpT=VsVdTrace1.T(1,tmpVds1idx).T;
tmpVs=VsVdTrace1.Vs(1,tmpVds1idx).Vs;
plot(tmpT,tmpVs,'-k')
title(['E_{K}=',num2str(Ek1),' mV V_{ds}^{out}= ',num2str(samplevds(1)), ' mV'])
ylabel('V_{s} (mV)')
ylim([-30,100]);

subplot(6,3,5)
%Ek2  Vds1
tmpVds2idx=find(VsVdTrace2.VdsOut(1,:)==samplevds(1));
tmpT=VsVdTrace2.T(1,tmpVds2idx).T;
tmpVs=VsVdTrace2.Vs(1,tmpVds2idx).Vs;
plot(tmpT,tmpVs,'-k')
title(['E_{K}=',num2str(Ek2),' mV V_{ds}^{out}= ',num2str(samplevds(1)), ' mV'])
ylabel('V_{s} (mV)')
ylim([-30,100]);

subplot(6,3,8)
%Ek1 vds2
tmpVds1idx=find(VsVdTrace1.VdsOut(1,:)==samplevds(2));
tmpT=VsVdTrace1.T(1,tmpVds1idx).T;
tmpVs=VsVdTrace1.Vs(1,tmpVds1idx).Vs;
plot(tmpT,tmpVs,'-k')
title(['E_{K}=',num2str(Ek1),' mV V_{ds}^{out}= ',num2str(samplevds(2)), ' mV'])
ylabel('V_{s} (mV)')
ylim([-30,100]);
%
subplot(6,3,11)
% Ek2 Vds2
tmpVds2idx=find(VsVdTrace2.VdsOut(1,:)==samplevds(2));
tmpT=VsVdTrace2.T(1,tmpVds2idx).T;
tmpVs=VsVdTrace2.Vs(1,tmpVds2idx).Vs;
plot(tmpT,tmpVs,'-k')
title(['E_{K}=',num2str(Ek2),' mV V_{ds}^{out}= ',num2str(samplevds(2)), ' mV'])
ylabel('V_{s} (mV)')
ylim([-30,100]);

%
subplot(6,3,14)
% Ek1 Vds3
tmpVds1idx=find(VsVdTrace1.VdsOut(1,:)==samplevds(3));
tmpT=VsVdTrace1.T(1,tmpVds1idx).T;
tmpVs=VsVdTrace1.Vs(1,tmpVds1idx).Vs;
plot(tmpT,tmpVs,'-k')
title(['E_{K}=',num2str(Ek1),' mV V_{ds}^{out}= ',num2str(samplevds(3)), ' mV'])
ylabel('V_{s} (mV)')
ylim([-30,100]);

subplot(6,3,17)
% Ek2 Vds3
tmpVds2idx=find(VsVdTrace2.VdsOut(1,:)==samplevds(3));
tmpT=VsVdTrace2.T(1,tmpVds2idx).T;
tmpVs=VsVdTrace2.Vs(1,tmpVds2idx).Vs;
plot(tmpT,tmpVs,'-k')
title(['E_{K}=',num2str(Ek2),' mV V_{ds}^{out}= ',num2str(samplevds(3)), ' mV'])
xlabel('Time (ms)')
ylabel('V_{s} (mV)')
ylim([-30,100]);

subplot(6,3,3)
%Ek1 Vds1
tmpVds1idx=find(VsVdTrace1.VdsOut(1,:)==samplevds(1));
tmpT=VsVdTrace1.T(1,tmpVds1idx).T;
tmpVs=VsVdTrace1.Vs(1,tmpVds1idx).Vs;
plot(tmpT,tmpVs,'-k')
title(['E_{K}=',num2str(Ek1),' mV V_{ds}^{out}= ',num2str(samplevds(1)), ' mV'])
ylabel('V_{s} (mV)')
ylim([-30,100]);

subplot(6,3,6)
%Ek2  Vds1
tmpVds2idx=find(VsVdTrace2.VdsOut(1,:)==samplevds(1));
tmpT=VsVdTrace2.T(1,tmpVds2idx).T;
tmpVs=VsVdTrace2.Vs(1,tmpVds2idx).Vs;
plot(tmpT,tmpVs,'-k')
title(['E_{K}=',num2str(Ek2),' mV V_{ds}^{out}= ',num2str(samplevds(1)), ' mV'])
ylabel('V_{s} (mV)')
ylim([-30,100]);

subplot(6,3,9)
%Ek1 vds2
tmpVds1idx=find(VsVdTrace1.VdsOut(1,:)==samplevds(2));
tmpT=VsVdTrace1.T(1,tmpVds1idx).T;
tmpVs=VsVdTrace1.Vs(1,tmpVds1idx).Vs;
plot(tmpT,tmpVs,'-k')
title(['E_{K}=',num2str(Ek1),' mV V_{ds}^{out}= ',num2str(samplevds(2)), ' mV'])
ylabel('V_{s} (mV)')
ylim([-30,100]);
%
subplot(6,3,12)
% Ek2 Vds2
tmpVds2idx=find(VsVdTrace2.VdsOut(1,:)==samplevds(2));
tmpT=VsVdTrace2.T(1,tmpVds2idx).T;
tmpVs=VsVdTrace2.Vs(1,tmpVds2idx).Vs;
plot(tmpT,tmpVs,'-k')
title(['E_{K}=',num2str(Ek2),' mV V_{ds}^{out}= ',num2str(samplevds(2)), ' mV'])
ylabel('V_{s} (mV)')
ylim([-30,100]);

%
subplot(6,3,15)
% Ek1 Vds3
tmpVds1idx=find(VsVdTrace1.VdsOut(1,:)==samplevds(3));
tmpT=VsVdTrace1.T(1,tmpVds1idx).T;
tmpVs=VsVdTrace1.Vs(1,tmpVds1idx).Vs;
plot(tmpT,tmpVs,'-k')
title(['E_{K}=',num2str(Ek1),' mV V_{ds}^{out}= ',num2str(samplevds(3)), ' mV'])
ylabel('V_{s} (mV)')
ylim([-30,100]);

subplot(6,3,18)
% Ek2 Vds3
tmpVds2idx=find(VsVdTrace2.VdsOut(1,:)==samplevds(3));
tmpT=VsVdTrace2.T(1,tmpVds2idx).T;
tmpVs=VsVdTrace2.Vs(1,tmpVds2idx).Vs;
plot(tmpT,tmpVs,'-k')
title(['E_{K}=',num2str(Ek2),' mV V_{ds}^{out}= ',num2str(samplevds(3)), ' mV'])
xlabel('Time (ms)')
ylabel('V_{s} (mV)')
ylim([-30,100]);

%Commented out below is just plots counting the number of spikes

% figure(); 
% subplot(2,2,1); 
% title(['E_{K}=',num2str(Ek1),' mV M= ',num2str(M),' \muA/(cm^{2} s)'])
% plot(VsVdTrace1.VdsOut(1,:),VsVdTrace1.NumSpks(1,:),'or');ylim([0,11]);
% subplot(2,2,3);
% title(['E_{K}=',num2str(Ek1),' mV M= ',num2str(M),' \muA/(cm^{2} s)'])
% semilogy(VsVdTrace1.VdsOut(1,:),VsVdTrace1.FirstToLastSpike(1,:));
% subplot(2,2,2); 
% title(['E_{K}=',num2str(Ek2),' mV M= ',num2str(M),' \muA/(cm^{2} s)'])
% plot(VsVdTrace2.VdsOut(1,:),VsVdTrace2.NumSpks(1,:),'or');ylim([0,11]); 
% subplot(2,2,4)
% title(['E_{K}=',num2str(Ek2),' mV M= ',num2str(M),' \muA/(cm^{2} s)'])
% subplot(2,2,4),semilogy(VsVdTrace2.VdsOut(1,:),VsVdTrace2.FirstToLastSpike(1,:))

% savefig(f1,'SpikeTypMpt3Ekm25m45.fig')
% save('VsTraceMpt3Ekm25.mat',VsVdTrace1);
% save('VsTraceMpt3Ekm45.mat',VsVdTrace2);
%% Loop over soecified polarizations and integrate the four user specified $E_{K}$ and $g_{KAHP}$ and plot soma potentials

