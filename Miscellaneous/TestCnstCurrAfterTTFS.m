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
Tend=8000;
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
M=0.3;
tmpuAmpsPermsecCm2=M*1e-3; %per msec
%NumtmpVdsOut=33;
%NumtmpVdsOut=15; %from -12.25 this goes to -16
%NumtmpVdsOut=75; %from -12.25 this goes to -16
NumVdsOut=12;
VdsTTFS1=zeros(NumVdsOut,2);
%StarttmpVdsOut=-4;
%StarttmpVdsOut=-12.25;
StarttmpVdsOut=-3;
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
samplevds=[0,-7,-13];
tafter=15000; %ms
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
   subsample=100;
   SingPRWithInjCurr1.T=[SingPRWithInjCurrbefore.T',SingPRWithInjCurrafter.T']';
   SingPRWithInjCurr1.YMultcol=[SingPRWithInjCurrbefore.YMultcol',SingPRWithInjCurrafter.YMultcol']';
   extr=1:subsample:size(SingPRWithInjCurr1.T,1);
   VsVdTrace1.T(i).T=SingPRWithInjCurr1.T(extr);
   VsVdTrace1.Vs(i).Vs=SingPRWithInjCurr1.YMultcol(extr,1);
   VsVdTrace1.Vd(i).Vd=SingPRWithInjCurr1.YMultcol(extr,2);
   VsVdTrace1.VdsOut(i)=tmpVdsOut;
   VsVdTrace1.TTFS(i)=SingPRWithInjCurrbefore.te(1)-SingPRWithInjCurrbefore.delay;
   VsVdTrace1.NumSpks(i)=size(SingPRWithInjCurrbefore.te,1)+size(SingPRWithInjCurrafter.te,1);  %DOes passing throug Vs on the way up before  mean that it will not be conte in after? Assume yes
   VsVdTrace1.FirstToLastSpike(i)=SingPRWithInjCurrafter.te(size(SingPRWithInjCurrafter.te,1))-SingPRWithInjCurrbefore.te(1);
   VsVdTrace1.M=M;
   VsVdTrace1.Ek=Ek1;
   % SingPRWithSyn.datetime
else
    stable1=0;
end
SomaInj=true;
if stable1 == 1
    if isempty(SingPRWithInjCurrbefore.te) 
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
end
samplevds=[-4,-9,-14];
figure()
subplot(3,1,1)
tmpVds1idx=find(VsVdTrace1.VdsOut(1,:)==samplevds(1));
tmpT=VsVdTrace1.T(1,tmpVds1idx).T;
tmpVs=VsVdTrace1.Vs(1,tmpVds1idx).Vs;
plot(tmpT,tmpVs,'-k')
subplot(3,1,2)
tmpVds1idx=find(VsVdTrace1.VdsOut(1,:)==samplevds(2));
tmpT=VsVdTrace1.T(1,tmpVds1idx).T;
tmpVs=VsVdTrace1.Vs(1,tmpVds1idx).Vs;
plot(tmpT,tmpVs,'-k')
subplot(3,1,3)
tmpVds1idx=find(VsVdTrace1.VdsOut(1,:)==samplevds(3));
tmpT=VsVdTrace1.T(1,tmpVds1idx).T;
tmpVs=VsVdTrace1.Vs(1,tmpVds1idx).Vs;
plot(tmpT,tmpVs,'-k')
