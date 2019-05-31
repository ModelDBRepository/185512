Ek45Mpt8=load('TTFSEk45Mpt8PlusIh5ghVdsTTFS.mat');
Ek25Mpt8=load('TTFSEk25Mpt8PlusIh5ghVdsTTFS.mat');
Ek45Mpt3=load('TTFSEk45Mpt3PlusIh5ghVdsTTFS.mat');
Ek25Mpt3=load('TTFSEk25Mpt3PlusIh5ghVdsTTFS.mat');

Lippertghvals=[0.0,0.03,0.035,0.04,0.06];
Lipperthvhalfvals=[-85,-81,-78,-75,-71];

symbolbag={'-ks','-k^','-ko',':ks',':k^',':ko','-kp'};

figure()
subplot(2,2,1)
for z=1:size(Lippertghvals,2)
    plot(Ek45Mpt8.ghVdsOutTTFS(z,:,2),Ek45Mpt8.ghVdsOutTTFS(z,:,3),symbolbag{z})
    hold on;
    strlegend{z}=['max g_{h}= ',num2str(Lippertghvals(1,z)),' mS/cm^{2} V_{1/2}= ',num2str(Lipperthvhalfvals(1,z)+aPRwH.DiffRefVoltHcurr),' mV'];
end
h=legend(strlegend)
xlabel('V_{ds}^{out} (mV)')
ylabel('TTFS (ms)')
title('M= 0.8 \muA/(cm^{2}s)  E_{K}= -45 mV')
subplot(2,2,2)
for z=1:size(Lippertghvals,2)
    plot(Ek25Mpt8.ghVdsOutTTFS(z,:,2),Ek25Mpt8.ghVdsOutTTFS(z,:,3),symbolbag{z})
    hold on;
end
%h=legend(strlegend)
xlabel('V_{ds}^{out} (mV)')
ylabel('TTFS (ms)')
title('M= 0.8 \muA/(cm^{2}s)  E_{K}= -25 mV')
subplot(2,2,3)
for z=1:size(Lippertghvals,2)
    plot(Ek45Mpt3.ghVdsOutTTFS(z,:,2),Ek45Mpt3.ghVdsOutTTFS(z,:,3),symbolbag{z})
    hold on;
end
%h=legend(strlegend)
xlabel('V_{ds}^{out} (mV)')
ylabel('TTFS (ms)')
title('M= 0.3 \muA/(cm^{2}s)  E_{K}= -45 mV')
subplot(2,2,4)
for z=1:size(Lippertghvals,2)
    plot(Ek25Mpt3.ghVdsOutTTFS(z,:,2),Ek25Mpt3.ghVdsOutTTFS(z,:,3),symbolbag{z})
    hold on;
end
%h=legend(strlegend)
xlabel('V_{ds}^{out} (mV)')
ylabel('TTFS (ms)')
title('M= 0.3 \muA/(cm^{2}s)  E_{K}= -25 mV')




