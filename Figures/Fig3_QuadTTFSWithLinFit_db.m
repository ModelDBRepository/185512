TTFSVdsIrk_Ismpt5=load('TTFSVdsIrIsm0pt5Km20nm40T12000Irpt05topt8persec');
VdsIrk_Ismpt5=load('VdsIrIsm0pt5Km20nm40T12000Irpt05topt8persec');
NumPot=11;
NumRamps=16;
Vk=zeros(NumPot-2,1);
m=zeros(NumRamps,1);
for pot=3:NumPot
Vk(pot)=-20-2.5*(pot-1); 
end
for k = 1:NumRamps
m(k)=0.05+.05*(k-1);
end
eps=0.025;
VdsMax=19;
VdsMin=-15;
% VdsMax=10;
% VdsMin=-12;
indexVdsMax=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMax-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMax+eps);
indexVdsMin=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMin-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMin+eps);
if indexVdsMax > indexVdsMin
    indexVds1=indexVdsMin;
    indexVds2=indexVdsMax;
else
     indexVds2=indexVdsMin;
    indexVds1=indexVdsMax;
end
mwant1=0.3;
indexmk1=find(m>mwant1-eps & m<mwant1+eps);
mwant2=0.8;
indexmk2=find(m>mwant2-eps & m<mwant2+eps);
Ekwant1=-25;
indexEk1=find(Vk>Ekwant1-eps & Vk<Ekwant1+eps);
Ekwant2=-45;
indexEk2=find(Vk>Ekwant2-eps & Vk<Ekwant2+eps);



%save(['TTFSIsmpt5Ir',num2str(VdsIrk_Ismpt5.VdsIrK(1,indexmk,5,2))],fig1)
fig2=figure()
subplot(2,2,1)
plot(VdsIrk_Ismpt5.VdsIrK(indexVds1:indexVds2,indexmk1,indexEk1,1),TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVds1:indexVds2,indexmk1,indexEk1),'k.')

title({['m = ', num2str(VdsIrk_Ismpt5.VdsIrK(1,indexmk1,indexEk1,2)*1e3),' Ek= ',num2str(VdsIrk_Ismpt5.VdsIrK(1,indexmk1,indexEk1,3))]})
eps=0.025;
s=0;
HitRsq=false;
 VdsMaxfit=5;
    VdsMinfit=0-(s-1)*1;
  
    indexVdsMaxfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMaxfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMaxfit+eps);
    indexVdsMinfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMinfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMinfit+eps);
    if indexVdsMaxfit > indexVdsMinfit
        indexVdsfit1=indexVdsMinfit;
        indexVdsfit2=indexVdsMaxfit;
    else
        indexVdsfit2=indexVdsMinfit;
        indexVdsfit1=indexVdsMaxfit;
    end

        [b,bint,r,rint,stats]=regress(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk1,1),[ones(size(VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk1))) VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk1)]);
beginRsq=1;
while beginRsq >=0.99 &VdsMinfit > VdsMaxfit-10%& s<10%& HitRsq==false
    s=s+1;
    VdsMaxfit=5;
    VdsMinfit=0-(s-1)*1;
    indexVdsMaxfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMaxfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMaxfit+eps);
    indexVdsMinfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMinfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMinfit+eps);
    if indexVdsMaxfit > indexVdsMinfit
        indexVdsfit1=indexVdsMinfit;
        indexVdsfit2=indexVdsMaxfit;
    else
        indexVdsfit2=indexVdsMinfit;
        indexVdsfit1=indexVdsMaxfit;
    end
           [b,bint,r,rint,stats]=regress(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk1,1),[ones(size(VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk1))) VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk1)]);
 
    beginRsq=stats(1);

end
    s=s-1;
    VdsMaxfit=5;
    VdsMinfit=0-(s-1)*1;
    indexVdsMaxfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMaxfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMaxfit+eps);
    indexVdsMinfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMinfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMinfit+eps);
    if indexVdsMaxfit > indexVdsMinfit
        indexVdsfit1=indexVdsMinfit;
        indexVdsfit2=indexVdsMaxfit;
    else
        indexVdsfit2=indexVdsMinfit;
        indexVdsfit1=indexVdsMaxfit;
    end
           [b,bint,r,rint,stats]=regress(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk1,1),[ones(size(VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk1))) VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk1)]);
  
    hline = refline([b(2) b(1)]);
    set(hline,'Color','r','LineStyle','-.','LineWidth',2)
    maxTTFS=max(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsMax:indexVdsMin,indexmk1,indexEk1,1))
    hold on;
    plot(VdsIrk_Ismpt5.VdsIrK(indexVdsfit2,indexmk1,indexEk1),TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit2,indexmk1,indexEk1,1),'LineStyle','None','MarkerSize',8,'Marker','*')
    ylim([0 maxTTFS*1.2])
      xlim([-20 12])
xlabel('V_{ds}^{out} (mv)')
ylabel('TTFS (ms)')
grid on;

subplot(2,2,4)
plot(VdsIrk_Ismpt5.VdsIrK(indexVds1:indexVds2,indexmk2,indexEk2,1),TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVds1:indexVds2,indexmk2,indexEk2),'k.')
title({['m = ', num2str(VdsIrk_Ismpt5.VdsIrK(1,indexmk2,indexEk2,2)*1e3),' Ek= ',num2str(VdsIrk_Ismpt5.VdsIrK(1,indexmk2,indexEk2,3))]})
eps=0.025;
s=0;
HitRsq=false;
 VdsMaxfit=5;
    VdsMinfit=0-(s-1)*1;
    indexVdsMaxfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMaxfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMaxfit+eps);
    indexVdsMinfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMinfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMinfit+eps);
    if indexVdsMaxfit > indexVdsMinfit
        indexVdsfit1=indexVdsMinfit;
        indexVdsfit2=indexVdsMaxfit;
    else
        indexVdsfit2=indexVdsMinfit;
        indexVdsfit1=indexVdsMaxfit;
    end
        [b,bint,r,rint,stats]=regress(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk2,1),[ones(size(VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk2))) VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk2)]);
beginRsq=1;
while beginRsq >=0.99 &VdsMinfit > VdsMaxfit-10%& s<10%& HitRsq==false
    s=s+1;
    VdsMaxfit=5;
    VdsMinfit=0-(s-1)*1;
    indexVdsMaxfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMaxfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMaxfit+eps);
    indexVdsMinfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMinfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMinfit+eps);
    if indexVdsMaxfit > indexVdsMinfit
        indexVdsfit1=indexVdsMinfit;
        indexVdsfit2=indexVdsMaxfit;
    else
        indexVdsfit2=indexVdsMinfit;
        indexVdsfit1=indexVdsMaxfit;
    end
           [b,bint,r,rint,stats]=regress(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk2,1),[ones(size(VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk2))) VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk2)]);
  
    beginRsq=stats(1);
end
    s=s-1;
    VdsMaxfit=5;
    VdsMinfit=0-(s-1)*1;
    indexVdsMaxfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMaxfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMaxfit+eps);
    indexVdsMinfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMinfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMinfit+eps);
    if indexVdsMaxfit > indexVdsMinfit
        indexVdsfit1=indexVdsMinfit;
        indexVdsfit2=indexVdsMaxfit;
    else
        indexVdsfit2=indexVdsMinfit;
        indexVdsfit1=indexVdsMaxfit;
    end
           [b,bint,r,rint,stats]=regress(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk2,1),[ones(size(VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk2))) VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk2)]);
  
    hline = refline([b(2) b(1)]);
    set(hline,'Color','r','LineStyle','-.','LineWidth',2)
    maxTTFS=max(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsMax:indexVdsMin,indexmk2,indexEk2,1))
    hold on;
    plot(VdsIrk_Ismpt5.VdsIrK(indexVdsfit2,indexmk2,indexEk2),TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit2,indexmk2,indexEk2,1),'LineStyle','None','MarkerSize',8,'Marker','*')
    ylim([0 maxTTFS*1.2])
      xlim([-20 12])
xlabel('V_{ds}^{out} (mv)')
ylabel('TTFS (ms)')
grid on;
subplot(2,2,3)
plot(VdsIrk_Ismpt5.VdsIrK(indexVds1:indexVds2,indexmk1,indexEk2,1),TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVds1:indexVds2,indexmk1,indexEk2),'k.')
title({['m = ', num2str(VdsIrk_Ismpt5.VdsIrK(1,indexmk1,indexEk2,2)*1e3),' Ek= ',num2str(VdsIrk_Ismpt5.VdsIrK(1,indexmk1,indexEk2,3))]})
eps=0.025;
s=0;
HitRsq=false;
 VdsMaxfit=5;
    VdsMinfit=0-(s-1)*1;
    indexVdsMaxfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMaxfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMaxfit+eps);
    indexVdsMinfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMinfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMinfit+eps);
    if indexVdsMaxfit > indexVdsMinfit
        indexVdsfit1=indexVdsMinfit;
        indexVdsfit2=indexVdsMaxfit;
    else
        indexVdsfit2=indexVdsMinfit;
        indexVdsfit1=indexVdsMaxfit;
    end

        [b,bint,r,rint,stats]=regress(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk2,1),[ones(size(VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk2))) VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk2)]);
beginRsq=1;
while beginRsq >=0.99 &VdsMinfit > VdsMaxfit-10%& s<10%& HitRsq==false
    s=s+1;
    VdsMaxfit=5;
    VdsMinfit=0-(s-1)*1;
    indexVdsMaxfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMaxfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMaxfit+eps);
    indexVdsMinfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMinfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMinfit+eps);
    if indexVdsMaxfit > indexVdsMinfit
        indexVdsfit1=indexVdsMinfit;
        indexVdsfit2=indexVdsMaxfit;
    else
        indexVdsfit2=indexVdsMinfit;
        indexVdsfit1=indexVdsMaxfit;
    end
           [b,bint,r,rint,stats]=regress(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk2,1),[ones(size(VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk2))) VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk2)]);
  
    beginRsq=stats(1);

end
    s=s-1;
    VdsMaxfit=5;
    VdsMinfit=0-(s-1)*1;
    indexVdsMaxfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMaxfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMaxfit+eps);
    indexVdsMinfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMinfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMinfit+eps);
    if indexVdsMaxfit > indexVdsMinfit
        indexVdsfit1=indexVdsMinfit;
        indexVdsfit2=indexVdsMaxfit;
    else
        indexVdsfit2=indexVdsMinfit;
        indexVdsfit1=indexVdsMaxfit;
    end
           [b,bint,r,rint,stats]=regress(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk2,1),[ones(size(VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk2))) VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk1,indexEk2)]);
  
    hline = refline([b(2) b(1)]);
    set(hline,'Color','r','LineStyle','-.','LineWidth',2)
    maxTTFS=max(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsMax:indexVdsMin,indexmk1,indexEk2,1))
    hold on;
    plot(VdsIrk_Ismpt5.VdsIrK(indexVdsfit2,indexmk1,indexEk2),TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit2,indexmk1,indexEk2,1),'LineStyle','None','MarkerSize',8,'Marker','*')
    ylim([0 maxTTFS*1.2])
      xlim([-20 12])
xlabel('V_{ds}^{out} (mv)')
ylabel('TTFS (ms)')
grid on;
subplot(2,2,2)
plot(VdsIrk_Ismpt5.VdsIrK(indexVds1:indexVds2,indexmk2,indexEk1,1),TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVds1:indexVds2,indexmk2,indexEk1),'k.')
title({['m = ', num2str(VdsIrk_Ismpt5.VdsIrK(1,indexmk2,indexEk1,2)*1e3),' Ek= ',num2str(VdsIrk_Ismpt5.VdsIrK(1,indexmk2,indexEk1,3))]})
eps=0.025;
s=0;
HitRsq=false;
 VdsMaxfit=5;
    VdsMinfit=0-(s-1)*1;
    indexVdsMaxfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMaxfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMaxfit+eps);
    indexVdsMinfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMinfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMinfit+eps);
    if indexVdsMaxfit > indexVdsMinfit
        indexVdsfit1=indexVdsMinfit;
        indexVdsfit2=indexVdsMaxfit;
    else
        indexVdsfit2=indexVdsMinfit;
        indexVdsfit1=indexVdsMaxfit;
    end
        [b,bint,r,rint,stats]=regress(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk1,1),[ones(size(VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk1))) VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk1)]);
beginRsq=1;
while beginRsq >=0.99 &VdsMinfit > VdsMaxfit-10%& s<10%& HitRsq==false
    s=s+1;
    VdsMaxfit=5;
    VdsMinfit=0-(s-1)*1;
    indexVdsMaxfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMaxfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMaxfit+eps);
    indexVdsMinfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMinfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMinfit+eps);
    if indexVdsMaxfit > indexVdsMinfit
        indexVdsfit1=indexVdsMinfit;
        indexVdsfit2=indexVdsMaxfit;
    else
        indexVdsfit2=indexVdsMinfit;
        indexVdsfit1=indexVdsMaxfit;
    end
           [b,bint,r,rint,stats]=regress(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk1,1),[ones(size(VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk1))) VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk1)]);
      beginRsq=stats(1);
end
    s=s-1;
    VdsMaxfit=5;
    VdsMinfit=0-(s-1)*1;
    indexVdsMaxfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMaxfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMaxfit+eps);
    indexVdsMinfit=find(VdsIrk_Ismpt5.VdsIrK(:,3,3,1)>VdsMinfit-eps & VdsIrk_Ismpt5.VdsIrK(:,3,3,1)<VdsMinfit+eps);
    if indexVdsMaxfit > indexVdsMinfit
        indexVdsfit1=indexVdsMinfit;
        indexVdsfit2=indexVdsMaxfit;
    else
        indexVdsfit2=indexVdsMinfit;
        indexVdsfit1=indexVdsMaxfit;
    end
           [b,bint,r,rint,stats]=regress(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk1,1),[ones(size(VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk1))) VdsIrk_Ismpt5.VdsIrK(indexVdsfit1:indexVdsfit2,indexmk2,indexEk1)]);
 
    hline = refline([b(2) b(1)]);
    set(hline,'Color','r','LineStyle','-.','LineWidth',2)
    maxTTFS=max(TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsMax:indexVdsMin,indexmk2,indexEk1,1))
    hold on;
    plot(VdsIrk_Ismpt5.VdsIrK(indexVdsfit2,indexmk1,indexEk2),TTFSVdsIrk_Ismpt5.TTFSVdsIrK(indexVdsfit2,indexmk2,indexEk1,1),'LineStyle','None','MarkerSize',8,'Marker','*')
    ylim([0 maxTTFS*1.2])
    xlim([-20 12])

xlabel('V_{ds}^{out} (mv)')
ylabel('TTFS (ms)')
grid on;
xsc=cell(1,1)
xsc=strcat('TTFSIsmpt5IrVarmEk','.fig');%,afilename(1,2),afilename(1,3),afilename(1,4));
saveas(fig2,xsc)

