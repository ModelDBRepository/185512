%% Just like Fig 17 with Fixed AMPA and "grid" plot of TTFS for g_KAHP and E_K BUT FOR PR Plus H-Current
% RIR October 3, 2015
data1=load('FixAMPApt3Ismpt5TTFSVaryEkgKAHPHiRespt075VdsTTFSEkKAHPsmallstep.mat');
data3=load('FixAMPApt3Ismpt5TTFSVaryEkgKAHPPRwHgh06Vm11VdsTTFSEkKAHP.mat');
%data3=load('FixAMPApt3Ismpt5TTFSVaryEkgKAHPPRwHgh03Vm21VdsTTFSEkKAHP.mat');
data2=load('FixAMPApt3Ismpt5TTFSVaryEkgKAHPPRwHgh03Vm21VdsTTFSEkKAHPs075.mat');
numKAHP1=size(data1.VdsTTFSgKAHPEk,1);
numEk1=size(data1.VdsTTFSgKAHPEk,2);
numVds1=size(data1.VdsTTFSgKAHPEk,3);
numdata1=size(data1.VdsTTFSgKAHPEk,4);
numKAHP2=size(data2.VdsTTFSgKAHPEk,1);
numEk2=size(data2.VdsTTFSgKAHPEk,2);
numVds2=size(data2.VdsTTFSgKAHPEk,3);
numdata2=size(data2.VdsTTFSgKAHPEk,4);
numKAHP3=size(data3.VdsTTFSgKAHPEk,1);
numEk3=size(data3.VdsTTFSgKAHPEk,2);
numVds3=size(data3.VdsTTFSgKAHPEk,3);
numdata3=size(data3.VdsTTFSgKAHPEk,4);

maxvds=-12; %mV Be Careful think twice about chainging have to worry about return spiking after maximum
vdsslices1=0:-0.075:maxvds;
vdsslices2=0:-0.25:maxvds;
vdsslices3=0:-0.25:maxvds;
ttfsatslices1= AMPAPolarPRSolnAnalyzerSlices_db(vdsslices1,data1.VdsTTFSgKAHPEk);
ttfsatslices2= AMPAPolarPRSolnAnalyzerSlices_db(vdsslices2,data2.VdsTTFSgKAHPEk);
ttfsatslices3= AMPAPolarPRSolnAnalyzerSlices_db(vdsslices3,data3.VdsTTFSgKAHPEk);
%idxslices=find(VdsTTFSgKAHPEk(1,1,:,1)==vdsslices);
slices1=size(ttfsatslices1,3);
slices2=size(ttfsatslices2,3);
slices3=size(ttfsatslices3,3);

Eks1=data1.VdsTTFSgKAHPEk(:,:,1,4)';
Eks1=repmat(Eks1,1,1,1);
Eks2=data2.VdsTTFSgKAHPEk(:,:,1,4)';
Eks2=repmat(Eks2,1,1,1);
Eks3=data3.VdsTTFSgKAHPEk(:,:,1,4)';
Eks3=repmat(Eks3,1,1,1);
%Eks=reshape(Eks,size(VdsTTFSgKAHPEk(:,1,1,1),1),size(VdsTTFSgKAHPEk(1,:,1,1),2),slices);
KAHPs1=data1.VdsTTFSgKAHPEk(:,:,1,3)';
KAHPs1=repmat(KAHPs1,1,1,1);
KAHPs2=data2.VdsTTFSgKAHPEk(:,:,1,3)';
KAHPs2=repmat(KAHPs2,1,1,1);
KAHPs3=data3.VdsTTFSgKAHPEk(:,:,1,3)';
KAHPs3=repmat(KAHPs3,1,1,1);
%=reshape(KAHPs,size(VdsTTFSgKAHPEk(:,1,1,1),1),size(VdsTTFSgKAHPEk(1,:,1,1),2),slices);

%figure()

% for i=1:slices
%    
%    pcolor(KAHPs,Eks,ttfsatslices(:,:,i))
%    caxis([0,100])
%  %  colormap(gray(10))
% %    map = colormap; % current colormap
%    s1=linspace(.1,.8,5);
%    map=[s1;s1;s1]';
%    map=brighten(map,0.5);
%    colormap(map)
%    colorbar
%    title(['V_{ds}^{out} =', num2str(vdsslices(i)), ' mV'])
%    xlabel('E_{K} (mV)')
%    ylabel('g_{KAHP} (ms/cm^{2})')
%    ylim([0,1.0]);
%    F(i)=getframe(gcf);
% end



f1=figure();

subplot(3,1,1)
pcolor(KAHPs1,Eks1,ttfsatslices1(:,:,slices1))
s1=linspace(.1,.8,5);
map=[s1;s1;s1]';
map=brighten(map,0.5);
colormap(map)
colorbar
title(['(b) Maximum TTFS for V_{ds}^{out} \in [', num2str(vdsslices1(slices1)), ', 0 mV]'])
xlabel('E_{K} (mV)')
ylabel('g_{KAHP} (ms/cm^{2})')
ylim([0,1.0]);

subplot(3,1,2)
pcolor(KAHPs2,Eks2,ttfsatslices2(:,:,slices2))
s1=linspace(.1,.8,5);
map=[s1;s1;s1]';
map=brighten(map,0.5);
colormap(map)
colorbar
title(['(b) Maximum TTFS for V_{ds}^{out} \in [', num2str(vdsslices2(slices2)), ', 0 mV]'])
xlabel('E_{K} (mV)')
ylabel('g_{KAHP} (ms/cm^{2})')
ylim([0,1.0]);

subplot(3,1,3)
pcolor(KAHPs3,Eks3,ttfsatslices3(:,:,slices3))
s1=linspace(.1,.8,5);
map=[s1;s1;s1]';
map=brighten(map,0.5);
colormap(map)
colorbar
title(['(b) Maximum TTFS for V_{ds}^{out} \in [', num2str(vdsslices3(slices3)), ', 0 mV]'])
xlabel('E_{K} (mV)')
ylabel('g_{KAHP} (ms/cm^{2})')
ylim([0,1.0]);


