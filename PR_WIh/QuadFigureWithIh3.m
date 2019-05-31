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
Tend=10000;
SomaInj=true; %Optionally the current ramp may be injected in the dendrite for gc=2.1 and p=0.5 there is not much qualitative  irrelevant as long as $M=0$
%SomaInj=false;
termaftevent=true; % Stops after Vs first rises up and passes through VsThresh
gAMPA=0.0;
aSyn.gAMPA=gAMPA;
gNMDA=0;
guessVsVd=[0,0];
aSyn.gNMDA=gNMDA;
M=0.8; %uA/(cm^2 s)
tmpuAmpsPermsecCm2=M*1e-3; %per msec
Ek=-25.0;
VsPreSyn.Amp=30;
VsPreSyn.t_initexc=delay;
VsPreSyn.t_spkdur=1.2; %msec
Isinj=-0.5;
%Ek=-38.56; % mV

%gh=0.03; % mS/cm^2 control Lippert 2009
gh=0.0;
h_Vhalf=-81; % mV control Lippert 2009
aPRwH=IniPRwH_db(Isinj,Ek,gh,h_Vhalf);
guessVsVd=[0,0];
ghss=[0.0,0.03];
Lippertghvals=[0.0,0.03,0.035,0.04,0.06];
Lipperthvhalfvals=[-85,-81,-78,-75,-71];

for z=1:size(Lippertghvals,2)
    gh=Lippertghvals(1,z);
    h_Vhalf=Lipperthvhalfvals(1,z); % mV control Lippert 2009
    aPRwH=IniPRwH_db(Isinj,Ek,gh,h_Vhalf);
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
ghVdsOutTTFS(z,i,1)=gh;
ghVdsOutTTFS(z,i,2)=VdsOut;
ghVdsOutTTFS(z,i,3)=SingPRPlusHWithSyn.te-delay;
end
end
figure()
symbolbag={'-ks','-k^','-ko',':ks',':k^',':ko','-kp'};
for z=1:size(Lippertghvals,2)
    plot(ghVdsOutTTFS(z,:,2),ghVdsOutTTFS(z,:,3),symbolbag{z})
    hold on;
end
h=legend(strlegend)
xlabel('V_{ds}^{out} (mV)')
ylabel('TTFS (ms)')
save('TTFSEk25Mpt8PlusIh5ghVdsTTFS.mat','ghVdsOutTTFS')