%% Automate curvature or sublinear and superlinear detector
%RIR 7/14/2015


%% What it does:
% This routine uses the a numerical second derivative of TTFS wrt
% $V_{ds}^{out}$ to determine if a profile is sublinear or superlinear. An
% estimate of the error in each TTFS computation is used to compute an
% error in the calculated second derivative. 

%% Approach:
% Starting at the highest (or least negative $V_{ds}^{out}$ we check at
% either the resolution of the input data (i.e. delta $V_{ds}^{out}$) or we
% could subsample to smooth out the computations. If the second derivateve
% exceeds the error tolerance before either some specified $V_{ds}^{out}$
% or before the $V_{ds}^{out}$ corresponding to the maximum TTFS

%% Notes
% For many of the chosen paramters there will be regions of polarization
% for which the model neuron fails to spike

%% Output
% Output will be a 2 dimensional matrix one for each of the varied
% paramters. The value will be 1 for superlinear and -1 sublinear



function ttfsatslices= AMPAPolarPRSolnAnalyzerSlices_db(vdsslices,VdsTTFSP1P2)

numvds=size(VdsTTFSP1P2(1,1,:,1),3);
deltavds=VdsTTFSP1P2(1,1,2,1)-VdsTTFSP1P2(1,1,1,1);
%stepvds=deltavds*subsamp;
NumEk=size(VdsTTFSP1P2(1,:,1,4),2);
NumKAHP=size(VdsTTFSP1P2(:,1,1,3));
%vdsslices=0:-0.25:-12;%linspace(-12,0,);%[-4,-6,-8,-10,-12];
for z=1:NumEk
for j=1:NumKAHP % 12 gKAHP
    vds=VdsTTFSP1P2(j,z,:,1);
    numvd=size(vds,3);
    vds=reshape(vds,numvd,1);
    ttfs=VdsTTFSP1P2(j,z,:,2);
    ttfs=reshape(ttfs,numvd,1);
    idxNaN=find(isnan(ttfs));
    VdsAtNaN=vds(idxNaN);
    if isempty(VdsAtNaN)
        idxminNaN=NaN;
    else
        idxminVdsNaN=VdsAtNaN(1);
        idxminNaN=idxNaN(1);
    end
    [maxttfs,idxmaxttfs]=max(ttfs);
    vdsatmax=vds(idxmaxttfs);
    maxvds=-20;
    idxVds15=find(vds<=maxvds,1,'first');
    idxmaxvds=min([idxmaxttfs,idxVds15,idxminNaN]);
    exidx=find(vds<-4,1,'first');
    for s=1:size(vdsslices,2)
    idxslices(s)=find(vds<vdsslices(s),1,'first');
    end
    ttfsatslices(z,j,:)=ttfs(idxslices);
    
    
end
end
