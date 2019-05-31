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
Tend=500;
%SomaInj=true; %Optionally the current ramp may be injected in the dendrite for gc=2.1 and p=0.5 there is not much qualitative  irrelevant as long as $M=0$
SomaInj=false;
termaftevent=true; % Stops after Vs first rises up and passes through VsThresh
gAMPA=0.3;
aSyn.gAMPA=gAMPA;
gNMDA=0;
guessVsVd=[0,0];
aSyn.gNMDA=gNMDA;
M=0; %uA/(cm^2 s)
tmpuAmpsPermsecCm2=M*1e-3; %per msec
%Ek=-45.0;
VsPreSyn.Amp=30;
VsPreSyn.t_initexc=delay;
VsPreSyn.t_spkdur=1.2; %msec
Isinj=-0.5;


% Eks=[-25,-25,-35,-35];
% gKAHPs=[0.2,0.9,0.1,0.4];

Eks=[-27.5,-27.5,-37.5,-37.5];
gKAHPs=[0.4,0.9,0.2,0.5];

Lippertghvals=[0.0,0.03,0.035,0.04,0.06];
Lipperthvhalfvals=[-85,-81,-78,-75,-71];

symbolbag={'-ks','-k^','-ko',':ks',':k^',':ko','-kp'};

%aPRwH=IniPRwH_db(Isinj,Ek,gh,h_Vhalf);

ghVdsOutTTFS=zeros(4,5,40,2);
tic()
for w=1:4
   

for z=1:size(Lippertghvals,2)
    gh=Lippertghvals(1,z);
    h_Vhalf=Lipperthvhalfvals(1,z); % mV control Lippert 2009
    Ek=Eks(1,w);
    aPRwH=IniPRwH_db(Isinj,Ek,gh,h_Vhalf);
    aPRwH.gKAHP=gKAHPs(1,w);
    strlegend{z}=['max g_{h}= ',num2str(Lippertghvals(1,z)),' mS/cm^{2} V_{1/2}= ',num2str(Lipperthvhalfvals(1,z)+aPRwH.DiffRefVoltHcurr),' mV'];
for i=1:40
VdsOut=5-(i-1)*0.5;
[numSSPRwH,diffProjFullEq,Jacob,eigJacob,nzeig] = NumerEquilPRwHcurr_db(aPRwH,guessVsVd,VdsOut);
aPRwH.SS.NumSS=numSSPRwH;
aPRwH.SS.eig = eig(Jacob);
if aPRwH.SS.eig < 0 %check if all eigenvalues less than zero
    stable=1;
    SingPRPlusHWithSyn = SingIntegODE23PRPlusHcurrWithSyn_db(aPRwH,tmpuAmpsPermsecCm2,VsPreSyn,aSyn,delay,Tend,VdsOut,VsThresh,SomaInj,termaftevent);
    %            [curr, pot, Spikes] = getcurrForSynInput_db(SingPRWithSyn,aSyn,tmpVdsOut,tmpuAmpsPermsecCm2,delay);
else
    stable=0;
end
z
i
%ghVdsOutTTFS(w,z,i,1)=gh;
ghVdsOutTTFS(w,z,i,1)=VdsOut;
if isempty(SingPRPlusHWithSyn.te)
    ghVdsOutTTFS(w,z,i,2)=NaN;
else
    ghVdsOutTTFS(w,z,i,2)=SingPRPlusHWithSyn.te-delay;
end
end
end
end
figure()
for i=1:size(ghVdsOutTTFS,1)
    for j=1:size(ghVdsOutTTFS,2)
        tmpVdsOut=squeeze(ghVdsOutTTFS(i,j,:,1));
        tmpTTFS=squeeze(ghVdsOutTTFS(i,j,:,2));
        plot(tmpVdsOut,tmpTTFS)
        hold on;
    end
    z
    toc
end

etime=toc()