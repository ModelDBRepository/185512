load('C:\Users\robert\My Documents\MATLAB\AMPATTFSp0pt5Ismpt3VdsTTFSgKAHPgAMPAForSplitpt25.mat')

%% VdsTTFSgKAHPAMPA size (10x8x140x4) (NumKAHP,NumAMPA,NumVdsOut,Num Dim), Simensions (VdsOut,TTFS,
tmpVdsTTFSgKAHPAMPA=VdsTTFSgKAHPAMPA;
gKAHPs=[0.0,0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6];
gAMPAs=[0.15,0.18,0.21,0.24,0.27,0.3,0.33,0.36];
Ek=-38.56;
tmpVdsTTFSgKAHPAMPA(:,:,:,3)=repmat(gKAHPs',1,8,140); 

%% Make 8 panels each panel 1 of 8 g_AMPA
% Each panel with constant g_AMPA and Ek=-38.56 mV will have 10 series for
% each of the g_KAHPs. X-axis VdsOut and y-axis TTFS which are
% Data(:,:,:,1) and Data(:,:,:,2) respectively
symbag={'-ko','-ks','-kv','-k^','-k*','-k+','-kv','-kd',':k','-k'};
figure()
for z=1:size(tmpVdsTTFSgKAHPAMPA,2)
subplot(4,2,z)
for i=1:size(tmpVdsTTFSgKAHPAMPA,1)
    tmpVds=tmpVdsTTFSgKAHPAMPA(i,z,:,1);
    tmpTTFS=tmpVdsTTFSgKAHPAMPA(i,z,:,2);
    tmpVds=reshape(tmpVds,size(tmpVds,3),1);
    tmpTTFS=reshape(tmpTTFS,size(tmpTTFS,3),1);
    plot(tmpVds,tmpTTFS,symbag{i})
    hold on;
end
title(['g_{AMPA}= ',num2str(gAMPAs(1,z)),' \muA/(cm^{2}s)'],'FontSize',8)%' E_{K}= ',num2str(Ek),' mV'])
xlim([-30,0])
ylim([0,100])
end