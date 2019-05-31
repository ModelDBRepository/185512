 function [ProjQuasiEqVsVd,quasiSSVsVdCahnscqi,fval,exitflag,output,Jacob] = ProjGateEquilPRwHcurrVsVd_db(aPRwH,Vs,Vd,VdsOut)
    % May 23, 2015
 %% What does it do?
 % It evaluates all the gating variables at their equilibrium values thus
 % reducing the equilibrium to a function of just Vs and Vd. It then calls
 % MATLAB's fsolve to solve for the equilibrium Vs and Vd given some starting Vs and
 % Vd
 
    hInf=GetGateInf_db(alphah_db(Vs,aPRwH.WRT),betah_db(Vs,aPRwH.WRT));
    nInf=GetGateInf_db(alphan_db(Vs,aPRwH.WRT),betan_db(Vs,aPRwH.WRT));
    sInf=GetGateInf_db(alphas_db(Vd,aPRwH.WRT),betas_db(Vd,aPRwH.WRT));
    cInf=GetGateInf_db(alphac_db(Vd,aPRwH.WRT),betac_db(Vd,aPRwH.WRT));
    ICad=aPRwH.gCa*(Vd-aPRwH.ECa)*sInf^2;
    Ca=-0.13*ICad/0.075;
    qInf=GetGateInf_db(alphaq_db(Ca),betaq_db);
    mInf=GetGateInf_db(alpham_db(Vs,aPRwH.WRT),betam_db(Vs,aPRwH.WRT));
    iinf=GateEquil_db(alpha_idB(Vd,aPRwH.h_Vhalf,aPRwH.HcurrLippertWRT),beta_idB(Vd,aPRwH.h_Vhalf,aPRwH.HcurrLippertWRT));
    Chid = min(Ca/250,1);
    VsVd(1)=Vs;
    VsVd(2)=Vd;
    %options=optimset('TolFun',1e-14,'Display','iter');
    options=optimset('TolFun',1e-14,'Display','none');
    [ProjQuasiEqVsVd,fval,exitflag,output,Jacob] =fsolve(@GetProjPRwHVsVd,VsVd,options);
    %[ProjQuasiEqVsVd,fval,exitflag,Jacob] =fsolve(@GetProjVsVd,VsVd);
    Vs=ProjQuasiEqVsVd(1);
    Vd=ProjQuasiEqVsVd(2);
    quasiSSVsVdCahnscqi(1)=ProjQuasiEqVsVd(1);
    quasiSSVsVdCahnscqi(2)=ProjQuasiEqVsVd(2);
    quasiSSVsVdCahnscqi(3)=(-0.13*aPRwH.gCa*(Vd-aPRwH.ECa)*sInf^2)/0.075;
    quasiSSVsVdCahnscqi(4)=GetGateInf_db(alphah_db(Vs,aPRwH.WRT),betah_db(Vs,aPRwH.WRT));
    quasiSSVsVdCahnscqi(5)=GetGateInf_db(alphan_db(Vs,aPRwH.WRT),betan_db(Vs,aPRwH.WRT));
    quasiSSVsVdCahnscqi(6)=GetGateInf_db(alphas_db(Vd,aPRwH.WRT),betas_db(Vd,aPRwH.WRT));
    quasiSSVsVdCahnscqi(7)=GetGateInf_db(alphac_db(Vd,aPRwH.WRT),betac_db(Vd,aPRwH.WRT));
    quasiSSVsVdCahnscqi(8)=GetGateInf_db(alphaq_db(Ca),betaq_db);
    quasiSSVsVdCahnscqi(9)=GateEquil_db(alpha_idB(Vd,aPRwH.h_Vhalf,aPRwH.HcurrLippertWRT),beta_idB(Vd,aPRwH.h_Vhalf,aPRwH.HcurrLippertWRT));;
    
   
    function ProjPRwHVsVd=GetProjPRwHVsVd(VsVd)
    
        Vs=VsVd(1);
        Vd=VsVd(2);
        ProjPRwHVsVd(1)=(1/aPRwH.Cm)*(-aPRwH.gL*Vs-aPRwH.gNa*mInf^2*hInf*(Vs-aPRwH.ENa)-aPRwH.gKDR*nInf*(Vs-aPRwH.Ek)...
         + (aPRwH.gc/aPRwH.p)*(Vd+VdsOut-Vs)+aPRwH.Isinj/aPRwH.p);
        ProjPRwHVsVd(2)=(1/aPRwH.Cm)*(-aPRwH.gL*Vd-ICad-aPRwH.gKAHP*qInf*(Vd-aPRwH.Ek)-aPRwH.gKC*cInf*Chid*(Vd-aPRwH.Ek)...
         - aPRwH.gh*iinf*(Vd-aPRwH.Eh)+(aPRwH.gc/(1-aPRwH.p))*(Vs-VdsOut-Vd)+aPRwH.Idinj/(1-aPRwH.p));
    end

end