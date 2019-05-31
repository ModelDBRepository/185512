%% Script to test finding equilibrium for PR wHcurr model
% RIR 10.1.2015

Isinj=-0.5;
Ek=-38.56; % mV
gh=0.03; % mS/cm^2 control Lippert 2009
h_Vhalf=-81; % mV control Lippert 2009
aPRwH=IniPRwH_db(Isinj,Ek,gh,h_Vhalf);
guessVsVd=[0,0];
VdsOut=0;
ghvals=[0.0,0.015,0.03,0.045,0.06,0.3,0.6];

%% Below are gh and h_Vhalf paired values baseon on Lippert 2009. Note these are the single compartment values. The have two compartments but thatis with the H-current in the soma and dendrite passive
% Lippertghvals=[0.0,0.03,0.035,0.04,0.06];
% Lipperthvhalfvals=[-85,-81,-78,-75,-71];
Lippertghvals=[0.0,0.03,0.04,0.06];
Lipperthvhalfvals=[-85,-81,-75,-71];
figure()
plot(Lippertghvals',Lipperthvhalfvals')
xlabel('g_{h} (mS/cm^{2})')
ylabel('V_{1/2} (mV)')
title('Paired values of max g_{h} and V_{1/2} for 5 regulated I_{h} currents. From Lippert 2009') 
VdsOutVsVdi=zeros(size(Lippertghvals,2),30,5);
for k=1:size(Lippertghvals,2)
    tmpgh=Lippertghvals(1,k);
    tmphvhalf=Lipperthvhalfvals(1,k);
    aPRwH=IniPRwH_db(Isinj,Ek,tmpgh,tmphvhalf);
    if k==1
        strlegend{k}='Polarized PR No H-Current';
    else
        strlegend{k}=['max g_{h}= ',num2str(Lippertghvals(1,k)),' mS/cm^{2} V_{1/2}= ',num2str(Lipperthvhalfvals(1,k)+aPRwH.DiffRefVoltHcurr),' mV'];
    end
    for i=1:30
    i
    tmpVdsOut=-25+i;
    [numSSPRwH,diffProjFullEq,Jacob,eigJacob,nzeig] = NumerEquilPRwHcurr_db(aPRwH,guessVsVd,tmpVdsOut);
    zx=9;
    if nzeig==0
        VdsOutVsVdi(k,i,1)=tmpgh;
        VdsOutVsVdi(k,i,2)=tmpVdsOut;
        VdsOutVsVdi(k,i,3)=numSSPRwH(1,1);
        VdsOutVsVdi(k,i,4)=numSSPRwH(1,2);
        VdsOutVsVdi(k,i,5)=numSSPRwH(1,9);
    else
        VdsOutVsVdi(k,i,1)=tmpgh;
        VdsOutVsVdi(k,i,2)=tmpVdsOut;
        VdsOutVsVdi(k,i,3)=NaN;
        VdsOutVsVdi(k,i,4)=NaN;
        VdsOutVsVdi(k,i,5)=NaN;
    end
end
end
symbolbag={'-ks','-k^','-ko',':ks',':k^',':ko','-kp'};

figure()
subplot(3,1,1)
for k=1:size(Lippertghvals,2)
plot(VdsOutVsVdi(k,:,2),VdsOutVsVdi(k,:,3),symbolbag{k})
hold on;
end
% plot(VdsOutVsVdi(2,:,2),VdsOutVsVdi(2,:,3),symbolbag{2})
% hold on;
% plot(VdsOutVsVdi(3,:,2),VdsOutVsVdi(3,:,3),symbolbag{3})
% hold on;
% plot(VdsOutVsVdi(4,:,2),VdsOutVsVdi(4,:,3),symbolbag{4})
% hold on;
% plot(VdsOutVsVdi(5,:,2),VdsOutVsVdi(5,:,3),symbolbag{5})
% hold on;
% plot(VdsOutVsVdi(6,:,2),VdsOutVsVdi(6,:,3),symbolbag{6})
% hold on;
% plot(VdsOutVsVdi(7,:,2),VdsOutVsVdi(7,:,3),symbolbag{7})
h=legend(strlegend)
%set(h,'Interpreter','latex','fontsize',24)


%% Add just PR for sanity check should match PRwHcurr with gh=0

% aPR=IniPR_db(Isinj,Ek);
% guessVsVd=[0,0];
% VdsOutVsVd=zeros(100,3);
% for i = 1:160
% 
% tmpVdsOut=-40+i*0.5;
% 
% [numSSPR,diffProjFullEq,Jacob,eigJacob,nzeig] = NumerEquilPR_db(aPR,guessVsVd,tmpVdsOut);
% if nzeig==0
%     VdsOutVsVd(i,1)=tmpVdsOut;
%     VdsOutVsVd(i,2)=numSSPR(1);
%     VdsOutVsVd(i,3)=numSSPR(2);
% else
%     VdsOutVsVd(i,1)=tmpVdsOut;
%     VdsOutVsVd(i,2)=NaN;
%     VdsOutVsVd(i,3)=NaN;
% end
% aPR.SS.NumSS=numSSPR;
% aPR.SS.eig = eig(Jacob);
% 
% end
% hold on
% plot(VdsOutVsVd(:,1),VdsOutVsVd(:,2),'rs')
% xlabel('V_{ds}^{out} (mV)')
% ylabel('V_{s} (mV)')
% figure()
% subplot(2,2,1)
% plot(VdsOutVsVdi(:,1),VdsOutVsVdi(:,2))
% title(['gh=',num2str(aPRwH.gh),' mS/cm^{2}'])
% xlabel('V_{DS}^{Out} (mV)')
% ylabel('V_{s} (mV)')
% subplot(2,2,2)
% plot(VdsOutVsVdi(:,1),VdsOutVsVdi(:,3))
% title(['gh=',num2str(aPRwH.gh),' mS/cm^{2}'])
% xlabel('V_{DS}^{Out} (mV)')
% ylabel('V_{d} (mV)')
% subplot(2,2,3)
% plot(VdsOutVsVdi(:,3),VdsOutVsVdi(:,4))
% title(['gh=',num2str(aPRwH.gh),' mS/cm^{2}'])
% xlabel('V_{d} (mV)')
% ylabel('i')
% subplot(2,2,4)
% plot(VdsOutVsVdi(:,2),VdsOutVsVdi(:,3))
% title(['gh=',num2str(aPRwH.gh),' mS/cm^{2}'])
% xlabel('V_{s} (mV)')
% ylabel('V_{d} (mV)')
subplot(3,1,2)
for k=1:size(Lippertghvals,2)
plot(VdsOutVsVdi(k,:,2),VdsOutVsVdi(k,:,4),symbolbag{k})
hold on;
end
% plot(VdsOutVsVdi(2,:,2),VdsOutVsVdi(2,:,3),symbolbag{2})
% hold on;
% plot(VdsOutVsVdi(3,:,2),VdsOutVsVdi(3,:,3),symbolbag{3})
% hold on;
% plot(VdsOutVsVdi(4,:,2),VdsOutVsVdi(4,:,3),symbolbag{4})
% hold on;
% plot(VdsOutVsVdi(5,:,2),VdsOutVsVdi(5,:,3),symbolbag{5})
% hold on;
% plot(VdsOutVsVdi(6,:,2),VdsOutVsVdi(6,:,3),symbolbag{6})
% hold on;
% plot(VdsOutVsVdi(7,:,2),VdsOutVsVdi(7,:,3),symbolbag{7})
h=legend(strlegend)
%set(h,'Interpreter','latex','fontsize',24)


% %% Add just PR for sanity check should match PRwHcurr with gh=0
% 
% aPR=IniPR_db(Isinj,Ek);
% guessVsVd=[0,0];
% VdsOutVsVd=zeros(100,3);
% for i = 1:160
% 
% tmpVdsOut=-40+i*0.5;
% 
% [numSSPR,diffProjFullEq,Jacob,eigJacob,nzeig] = NumerEquilPR_db(aPR,guessVsVd,tmpVdsOut);
% 
% if nzeig==0
%     VdsOutVsVd(i,1)=tmpVdsOut;
%     VdsOutVsVd(i,2)=numSSPR(1);
%     VdsOutVsVd(i,3)=numSSPR(2);
% else
%     VdsOutVsVd(i,1)=tmpVdsOut;
%     VdsOutVsVd(i,2)=NaN;
%     VdsOutVsVd(i,3)=NaN;
% end
% 
% aPR.SS.NumSS=numSSPR;
% aPR.SS.eig = eig(Jacob);
% 
% end
% hold on
% plot(VdsOutVsVd(:,1),VdsOutVsVd(:,3),'rs')
% xlabel('V_{ds}^{out} (mV)')
% ylabel('V_{d} (mV)')

subplot(3,1,3)
for k=1:size(Lippertghvals,2)
plot(VdsOutVsVdi(k,:,2),VdsOutVsVdi(k,:,5),symbolbag{k})
hold on;
end
xlabel('V_{ds}^{out} (mV)')
ylabel('h activation gating variable i')