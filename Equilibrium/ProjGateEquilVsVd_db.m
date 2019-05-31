 function [ProjQuasiEqVsVd,quasiSSVsVdCahnscq,fval,exitflag,output,Jacob] = ProjGateEquilVsVd_db(aPR,Vs,Vd,VdsOut)
    % May 23, 2015
 %% What does it do?
 % It evaluates all the gating variables at their equilibrium values thus
 % reducing the equilibrium to a function of just Vs and Vd. It then calls
 % MATLAB's fsolve to solve for the equilibrium Vs and Vd given some starting Vs and
 % Vd
 
    hInf=GetGateInf_db(alphah_db(Vs,aPR.WRT),betah_db(Vs,aPR.WRT));
    nInf=GetGateInf_db(alphan_db(Vs,aPR.WRT),betan_db(Vs,aPR.WRT));
    sInf=GetGateInf_db(alphas_db(Vd,aPR.WRT),betas_db(Vd,aPR.WRT));
    cInf=GetGateInf_db(alphac_db(Vd,aPR.WRT),betac_db(Vd,aPR.WRT));
    ICad=aPR.gCa*(Vd-aPR.ECa)*sInf^2;
    Ca=-0.13*ICad/0.075;
    qInf=GetGateInf_db(alphaq_db(Ca),betaq_db);
    mInf=GetGateInf_db(alpham_db(Vs,aPR.WRT),betam_db(Vs,aPR.WRT));
    Chid = min(Ca/250,1);
    VsVd(1)=Vs;
    VsVd(2)=Vd;
    %options=optimset('TolFun',1e-14,'Display','iter');
    options=optimset('TolFun',1e-14,'Display','none');
    [ProjQuasiEqVsVd,fval,exitflag,output,Jacob] =fsolve(@GetProjVsVd,VsVd,options);
    %[ProjQuasiEqVsVd,fval,exitflag,Jacob] =fsolve(@GetProjVsVd,VsVd);
    Vs=ProjQuasiEqVsVd(1);
    Vd=ProjQuasiEqVsVd(2);
    quasiSSVsVdCahnscq(1)=ProjQuasiEqVsVd(1);
    quasiSSVsVdCahnscq(2)=ProjQuasiEqVsVd(2);
    quasiSSVsVdCahnscq(3)=(-0.13*aPR.gCa*(Vd-aPR.ECa)*sInf^2)/0.075;
    quasiSSVsVdCahnscq(4)=GetGateInf_db(alphah_db(Vs,aPR.WRT),betah_db(Vs,aPR.WRT));
    quasiSSVsVdCahnscq(5)=GetGateInf_db(alphan_db(Vs,aPR.WRT),betan_db(Vs,aPR.WRT));
    quasiSSVsVdCahnscq(6)=GetGateInf_db(alphas_db(Vd,aPR.WRT),betas_db(Vd,aPR.WRT));
    quasiSSVsVdCahnscq(7)=GetGateInf_db(alphac_db(Vd,aPR.WRT),betac_db(Vd,aPR.WRT));
    quasiSSVsVdCahnscq(8)=GetGateInf_db(alphaq_db(Ca),betaq_db);
    
   
    function ProjVsVd=GetProjVsVd(VsVd)
    
        Vs=VsVd(1);
        Vd=VsVd(2);
        ProjVsVd(1)=(1/aPR.Cm)*(-aPR.gL*Vs-aPR.gNa*mInf^2*hInf*(Vs-aPR.ENa)-aPR.gKDR*nInf*(Vs-aPR.Ek)...
         + (aPR.gc/aPR.p)*(Vd+VdsOut-Vs)+aPR.Isinj/aPR.p);
        ProjVsVd(2)=(1/aPR.Cm)*(-aPR.gL*Vd-ICad-aPR.gKAHP*qInf*(Vd-aPR.Ek)-aPR.gKC*cInf*Chid*(Vd-aPR.Ek)...
         + (aPR.gc/(1-aPR.p))*(Vs-VdsOut-Vd)+aPR.Idinj/(1-aPR.p));
    end

end