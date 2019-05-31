%% Just like Fig 17 with Fixed AMPA and "grid" plot of TTFS for g_KAHP and E_K BUT FOR PR Plus H-Current
% RIR October 3, 2015
%data=load('FixAMPApt3Ismpt5TTFSVaryEkgKAHPPRwHgh06Vm11VdsTTFSEkKAHP.mat');
data=load('FixAMPApt3Ismpt5TTFSVaryEkgKAHPPRwHgh03Vm21VdsTTFSEkKAHP.mat');
numKAHP=size(data.VdsTTFSgKAHPEk,1);
numEk=size(data.VdsTTFSgKAHPEk,2);
numVds=size(data.VdsTTFSgKAHPEk,3);
numdata=size(data.VdsTTFSgKAHPEk,4);
figure()
maxvds=-12; %mV Be Careful think twice about chainging have to worry about return spiking after maximum
vdsslices=0:-0.25:maxvds;
ttfsatslices= AMPAPolarPRSolnAnalyzerSlices_db(vdsslices,data.VdsTTFSgKAHPEk);
%idxslices=find(VdsTTFSgKAHPEk(1,1,:,1)==vdsslices);
slices=size(ttfsatslices,3);
F(slices) = struct('cdata',[],'colormap',[]);

Eks=data.VdsTTFSgKAHPEk(:,:,1,4)';
Eks=repmat(Eks,1,1,1);
%Eks=reshape(Eks,size(VdsTTFSgKAHPEk(:,1,1,1),1),size(VdsTTFSgKAHPEk(1,:,1,1),2),slices);
KAHPs=data.VdsTTFSgKAHPEk(:,:,1,3)';
KAHPs=repmat(KAHPs,1,1,1);
%=reshape(KAHPs,size(VdsTTFSgKAHPEk(:,1,1,1),1),size(VdsTTFSgKAHPEk(1,:,1,1),2),slices);
for i=1:slices
   
   pcolor(KAHPs,Eks,ttfsatslices(:,:,i))
   caxis([0,100])
 %  colormap(gray(10))
%    map = colormap; % current colormap
   s1=linspace(.1,.8,5);
   map=[s1;s1;s1]';
   map=brighten(map,0.5);
   colormap(map)
   colorbar
   title(['V_{ds}^{out} =', num2str(vdsslices(i)), ' mV'])
   xlabel('E_{K} (mV)')
   ylabel('g_{KAHP} (ms/cm^{2})')
   ylim([0,1.0]);
   F(i)=getframe(gcf);
end

f5=figure
%movie(f1,F,1,10)
movie(f5,F,1,10)
implay(F)


f2=figure();

subplot(2,1,2)
pcolor(KAHPs,Eks,ttfsatslices(:,:,slices))
s1=linspace(.1,.8,5);
map=[s1;s1;s1]';
map=brighten(map,0.5);
colormap(map)
colorbar
title(['(b) Maximum TTFS for V_{ds}^{out} \in [', num2str(vdsslices(i)), ', 0 mV]'])
xlabel('E_{K} (mV)')
ylabel('g_{KAHP} (ms/cm^{2})')
ylim([0,1.0]);

subplot(2,1,1)
idxmaxvds=find(data.VdsTTFSgKAHPEk(2,4,:,1)<maxvds,1,'first');
i=[1,3,4,7,6,11];%KAHP
j=[1,5,3,5,10,5];%Ek
ax2=gca();
sets=[-37.5,0.8;-35.0,0.5;-32.5,0.9;-32.5,1];
v=3;
symbag={'-sk','-+k','-ok','-dk'};
for z=1:size(sets,1)
    tmpEk=sets(z,1);
    tmpKAHP=sets(z,2);
    idxEk=find(data.VdsTTFSgKAHPEk(1,:,1,3)==tmpEk,1,'first');
    idxKAHP=find(data.VdsTTFSgKAHPEk(:,1,1,4)==tmpKAHP,1,'first');
    plot(squeeze(data.VdsTTFSgKAHPEk(idxKAHP,idxEk,1:6:idxmaxvds,1)),squeeze(data.VdsTTFSgKAHPEk(idxKAHP,idxEk,1:6:idxmaxvds,2)),symbag{z})
    top=0.4-.05*z;
%     tmpEk=VdsTTFSgKAHPEk(i(z),j(z),1,3);
%     tmpKAHP=VdsTTFSgKAHPEk(i(z),j(z),1,4);
leg1{z}=['E_{K}= ',num2str(tmpEk),' g_{KAHP}= ', ...
        num2str(tmpKAHP),' ms/cm^{2}']; % TTFS ',num2str(max(ttfsatslices(j(z),i(z),:))),' ms'];
%     tb(z)=annotation('textbox', [0.2 top 0.025 0.05], 'String', ['E_{K}= ',num2str(tmpEk),' g_{KAHP}= ', ...
%         num2str(tmpKAHP),' ms/cm^{2} TTFS ',num2str(max(ttfsatslices(j(z),i(z),:))),' ms']); %lower left (0,0) upper right (1,1)      
    hold on;
end

legend(leg1{1},leg1{2},leg1{3},leg1{4})
title('(a) PR with control-state of H-current.')
xlabel('V_{ds}^{out} (mV)')
ylabel('TTFS (ms)')
