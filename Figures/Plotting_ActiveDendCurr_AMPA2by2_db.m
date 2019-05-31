%% RIR June 7, 2015 
%% What it does.
% Integrates the Polarized PR model With *AMPA* Input for
% two sets of user defined $E_{k}$ and $g_{AMPA}$. It thenplots the total
% active dendrite currents side by side.
%% Purpose:
% To illustrate if the the same qualitative split in active dendrite
% currents that accompanies sub and super linear profiles for the ramp
% injection occur for synaptic AMPA input.

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
Tend=1000;
SomaInj=true; %Optionally the current ramp may be injected in the dendrite for gc=2.1 and p=0.5 there is not much qualitative  irrelevant as long as $M=0$
%SomaInj=false;
termaftevent=true; % Stops after Vs first rises up and passes through VsThresh
%% Initialize PR neuron

Idinj=0;
Isinj=-0.5; %uA/cm^2

Eks=[-25,-45]; %mv
%vds=[-6,-7,-8,-9,-10];
vds=[-6,-8,-10,-12];
M=0; %uA/(cm^2 s)
%aSyn.gAMPA=0.05;
aSyn.gNMDA=0;
tmpuAmpsPermsecCm2=M*1e-3; %per msec
VsPreSyn.Amp=30;
VsPreSyn.t_initexc=delay;
VsPreSyn.t_spkdur=1.2; %msec
% NumVdsOut=45;
% StartVdsOut=0;
% StepVdsOut=-0.25;
gAMPA=0.3;
endx=50;
aSyn.gAMPA=gAMPA;
gNMDA=0;
guessVsVd=[0,0];
aSyn.gNMDA=gNMDA;
set(gca,'LineStyleOrder',{'-.o','-.s','-.+','-.^','-o','-s','-+','-^','-*'})
lineorderstr={'-ok','-sk','-vk','-^k',':ok',':sk',':vk',':^k',':*k'};
Facecolor={'black','none','black','none','none','black','none','none','black'};
f1=figure()
%sets=[-38,0.6;-28,0.5];
%sets=[-35,.85;-35,0.2;-38,0.6;-28,0.5];
%sets=[-38,0.6;-35,.85;-35,0.2;-28,0.5];
%sets=[-25.0,1.2;-25.0,0.95;-25.0,0.85;-25.0,0.7];
%sets=[-30,0.9;-30,0.75;-30,0.65;-30,0.5];
sets=[-25.0,0.2;-35,0.1;-25.0,0.9;-35,0.4];
%setlabel={'Ek25_a','Ek25_b','Ek25_c','Ek25_d'};
%setlabel={'Ek30_a','Ek30_b','Ek30_c','Ek30_d'};
setlabel={'Ek25_bif','Ek30_bif','Ek35_bif','Ek40_bif'};
%setlabel={'B','C','A','D'};
for i=1:size(sets,1)
    subplot(2,2,i)
    %Ek=Eks(1,i);
    Ek=sets(i,1);
    gKAHP=sets(i,2);
    aPR=IniPR_db(Isinj,Ek);
    aPR.gKAHP=gKAHP;
    i
    for j=1:size(vds,2)
        tmpVdsOut=vds(1,j);
        [numSSPR,diffProjFullEq,Jacob,eigJacob,nzeig] = NumerEquilPR_db(aPR,guessVsVd,tmpVdsOut);
        
        aPR.SS.NumSS=numSSPR;
        aPR.SS.eig = eig(Jacob);
        if aPR.SS.eig < 0 %check if all eigenvalues less than zero
            stable=1;
            SingPRWithSyn = SingIntegODE23PRWithSyn_db(aPR,tmpuAmpsPermsecCm2,VsPreSyn,aSyn,delay,Tend,tmpVdsOut,VsThresh,SomaInj,termaftevent)  
            [curr, pot, Spikes] = getcurrForSynInput_db(SingPRWithSyn,aSyn,tmpVdsOut,tmpuAmpsPermsecCm2,delay);
        else
            stable=0;
        end
        j
      start=1;
    endplot=floor(SingPRWithSyn.idxteVs*1);
%    sub=5000;
    numTs=size(curr.T,1);
    sub=100;
    subsampleT=start+1:sub:endplot;
    subsampleCurr=start:sub:endplot-1;
    %subsample=start:sub:endplot;
%     plot(curr.T(subsample,1),curr.TotalActiveDend(subsample,1),lineorderstr{j},'MarkerFaceColor',Facecolor{j},'LineWidth',1,'MarkerSize',5)  
%plot(curr.T(2:numTs,1),curr.TotalActiveDend(1:numTs-1,1),lineorderstr{j},'MarkerFaceColor',Facecolor{j},'LineWidth',1,'MarkerSize',2)   
%plot(curr.T(2:numTs,1),curr.TotalActiveDend(1:numTs-1,1),lineorderstr{j},'LineWidth',1,'MarkerSize',2)   
plot(curr.T(subsampleT,1)-delay,curr.TotalActiveDend(subsampleCurr,1),lineorderstr{j},'LineWidth',1,'MarkerSize',8,'MarkerFaceColor',Facecolor{j})   
    hold on;  
     legstr{j}=num2str(vds(j));   
    end
    ylabel('active dendrite currents \muA/cm^{2}')
ylim([-0.2,10])
xlim([-5,endx])
xlabel('Time (ms)')
grid on
title(['E_{k}=',num2str(Ek),' mv  g_{KAHP}= ',num2str(gKAHP),' mS/cm^{2}.  Parameter Set: ',num2str(setlabel{i})],'FontSize',10,'FontWeight','normal')
h=gca;
set(h,'FontSize',10)
legend(legstr)
end

