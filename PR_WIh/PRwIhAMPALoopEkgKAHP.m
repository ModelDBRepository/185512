%% Loops through Ek and gKAHP at a fixed gAMPA for a range of $V_{ds}^{out}$
%% RIR 7/12/2015 
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
%Tend=800;
Tend=500; %m/s much quicker than ramp injected currents at least for the M we had been using. Actually most are under 100 ms
SomaInj=false; %Optionally the current ramp may be injected in the dendrite for gc=2.1 and p=0.5 there is not much qualitative difference
                % This does not matter since for this syanptic stimulus we
                % will set M=0;
termaftevent=true; % Stops after Vs first rises up and passes through VsThresh
%% Initialize PR neuron

Idinj=0; % Nothing extra is injected in the dendrite besides the synaptic current
%Isinj=-0.3; %uA/cm^2This should be selected to match known spike intervals times with injected current. -0.5 is the one used in PR and Grahm, Cutridis and Hunter in Associative model 
Isinj=-0.5;
% chapter in Hippocampal Microcircuits use -0.5 as. The origins of -0.3 is somehat of a mistake I used that when I also had FSI going. Would like to run it again sometime with Is=-0.5fits several other computational models using synaptic input and the PR model (Need References) 
%Isinj=0;
%Isinj=-1.0;
Ek=-15; %mv
Lippertghvals=[0.0,0.03,0.035,0.04,0.06];
Lipperthvhalfvals=[-85,-81,-78,-75,-71];
% gh=Lippertghvals(1,5);
% h_Vhalf=Lipperthvhalfvals(1,5);
gh=Lippertghvals(1,2);
h_Vhalf=Lipperthvhalfvals(1,2);
aPRwH=IniPR_db(Isinj,Ek);
aPRwH=IniPRwH_db(Isinj,Ek,gh,h_Vhalf);
M=0; %uA/(cm^2 s)
aSyn.gAMPA=0.3;
aSyn.gNMDA=0;
tmpuAmpsPermsecCm2=M*1e-3; %per msec
VsPreSyn.Amp=30;
VsPreSyn.t_initexc=delay;
VsPreSyn.t_spkdur=1.2; %msec
%NumVdsOut=33;
%NumVdsOut=15; %from -12.25 this goes to -16
%NumVdsOut=75; %from -12.25 this goes to -16
NumVdsOut=267;
%NumVdsOut=64;
%StartVdsOut=-4;
%StartVdsOut=-12.25;
StartVdsOut=0;
%StartVdsOut=-16.25;
%StepVdsOut=-0.25;
StepVdsOut=-0.075;
% MinVdsOyt=strtVdsOut+NumVdsOut*StepVdsOut 
MinVdsOut=StartVdsOut+NumVdsOut*StepVdsOut;
% error VdsTTFSEkgAMPA=zeros(NumVdsOut,2); found May 29, 2015
%% Integration Parameters
RelTol=1e-4;
AbsTol=1e-6;
MaxStep=1e-1;
IntegParam=[RelTol,AbsTol,MaxStep];

%% We will loop over VdsOut 
% first calculating the numerical equilibrium (If it exists). 
% If no equilibrium is found than we do not integrate

guessVsVd(1)=-10; % To initiate finding the equilibrium. -10, -10 is as good a place as any not real sensitive
guessVsVd(2)=-10;
%Eks=[-15,-20,-25,-30,-35,-40,-45,-50];
%Eks=[-25,-30,-35,-40,-45];
Eks=[-15,-17.5,-20,-22.5,-25,-27.5,-30,-32.5,-35,-37.5,-40];
%gAMPAs=[0.15,0.18,0.21,0.24,0.27,0.3,0.33,0.36];
%gAMPAs=0.27;
gAMPAs=0.30;
%gKAHPs=[0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6];
%gKAHPs=[0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6];
%gKAHPs=[0.0,0.05,0.1,0.15,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7];
gKAHPs=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2];
numgAMPAs=size(gAMPAs,2);
numgKAHPs=size(gKAHPs,2);
numEks=size(Eks,2);
VdsTTFSgKAHPEk=zeros(numgKAHPs,numEks,NumVdsOut,4);

aSyn.gAMPA=gAMPAs;


for p=1:numgKAHPs
    aPRwH.gKAHP=gKAHPs(p);
for g=1:numEks
    aPRwH.Ek=Eks(g);
    tmpEk=aPRwH.Ek;
    idx=(p-1)*numEks+g;
    percdone=idx/(numgKAHPs*numEks)
    display(['% done =',num2str(percdone)])
for i = 1:NumVdsOut

tmpVdsOut = StartVdsOut+(StepVdsOut*i);

[numSSPRwH,diffProjFullEq,Jacob,eigJacob,nzeig] = NumerEquilPRwHcurr_db(aPRwH,guessVsVd,tmpVdsOut);
aPRwH.SS.NumSS=numSSPRwH;
aPRwH.SS.eig = eig(Jacob);
if aPRwH.SS.eig < 0 %check if all eigenvalues less than zero
    stable=1;
    SingPRPlusHWithSyn = SingIntegODE23PRPlusHcurrWithSyn_db(aPRwH,tmpuAmpsPermsecCm2,VsPreSyn,aSyn,delay,Tend,tmpVdsOut,VsThresh,SomaInj,termaftevent);
    %            [curr, pot, Spikes] = getcurrForSynInput_db(SingPRWithSyn,aSyn,tmpVdsOut,tmpuAmpsPermsecCm2,delay);
else
    stable=0;
end
SomaInj=true;
if stable == 1
    if isempty(SingPRPlusHWithSyn.te) 
        VdsTTFSgKAHPEk(p,g,i,1)=tmpVdsOut;
        VdsTTFSgKAHPEk(p,g,i,2)=NaN; %NaN didn't spike
        VdsTTFSgKAHPEk(p,g,i,3)=tmpEk;
        VdsTTFSgKAHPEk(p,g,i,4)=aPRwH.gKAHP;
    else
        VdsTTFSgKAHPEk(p,g,i,1)=tmpVdsOut;
        VdsTTFSgKAHPEk(p,g,i,2)=SingPRPlusHWithSyn.te-SingPRPlusHWithSyn.delay; 
        VdsTTFSgKAHPEk(p,g,i,3)=tmpEk;
        VdsTTFSgKAHPEk(p,g,i,4)=aPRwH.gKAHP;
    end
else
    VdsTTFSgKAHPEk(p,g,i,1)=tmpVdsOut; 
    VdsTTFSgKAHPEk(p,g,i,2)=-1; %unstable
    VdsTTFSgKAHPEk(p,g,i,3)=tmpEk;
    VdsTTFSgKAHPEk(p,g,i,4)=aPRwH.gKAHP;
end
end
end
end
toc
save FixAMPApt3Ismpt5TTFSVaryEkgKAHPPRwHgh03Vm21VdsTTFSEkKAHPs075.mat VdsTTFSgKAHPEk
