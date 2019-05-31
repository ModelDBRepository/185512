  function [numSSPR,diffProjFullEq,fval,exitflag,output,Jacob,eigJacob,nzeig] = GetFzeroSSPR_db(aPR,VdsOut,quasiSS)

  %% What does it do?
  % This serves as another check for the steady state by using fsolve on
  % the full 8-dimensional problem. It uses as a starting point the
  % equilibrium state found by projecting the problem to the Vs Vd plane in
  % ProjGateEquilVsVd_db. The differenxe between the two should be minimal
  % One advantage is we get the full 8-d linearization with eigenvalues.
  
    %options=optimset('TolFun',1e-14,'Display','iter','MaxIter',1000,'MaxFunEvals',100000);
    options=optimset('TolFun',1e-14,'MaxIter',1000,'MaxFunEvals',100000,'Display','none');
    [numSSPR,fval,exitflag,output,Jacob] = fsolve(@PR94rhs,quasiSS,options);
    eigJacob=eig(Jacob);
    nzeig=0;
    for i=1:size(eigJacob,1)
       if eigJacob(i) > 0 
%           display('Warning one of the eigenvalues of the Jacobian is positive')
           nzeig=nzeig+1;
       end
    end
    diffProjFullEq=numSSPR-quasiSS;
function dY = PR94rhs(Y)
  
dY = zeros(8,1);    % a column vector
dY(1) = (1/aPR.Cm)*(-aPR.gL*(Y(1)-aPR.EL)-aPR.gNa*MInfPR94(Y(1),aPR.WRT)*Y(4)*(Y(1)-aPR.ENa)-aPR.gKDR*Y(5)*(Y(1)-aPR.Ek) ...
    +(aPR.gc/aPR.p)*(Y(2)+VdsOut-Y(1))+aPR.Isinj/aPR.p);
dY(2) = (1/aPR.Cm)*(-aPR.gL*(Y(2)-aPR.EL)-aPR.gCa*Y(6)^2*(Y(2)-aPR.ECa)-aPR.gKAHP*Y(8)*(Y(2)-aPR.Ek)-aPR.gKC*Y(7)*Chi(Y(3))*(Y(2)-aPR.Ek)...
        + aPR.gc/(1-aPR.p)*(Y(1)-VdsOut-Y(2))+aPR.Idinj/(1-aPR.p));
dY(3) = -0.13*aPR.gCa*Y(6)^2*(Y(2)-aPR.ECa)-0.075*Y(3);
dY(4) = (GateEquil_db(alphah_db(Y(1),aPR.WRT),betah_db(Y(1),aPR.WRT))-Y(4))/GateTimeCnst_db(alphah_db(Y(1),aPR.WRT),betah_db(Y(1),aPR.WRT));
dY(5) = (GateEquil_db(alphan_db(Y(1),aPR.WRT),betan_db(Y(1),aPR.WRT))-Y(5))/GateTimeCnst_db(alphan_db(Y(1),aPR.WRT),betan_db(Y(1),aPR.WRT));
dY(6) = (GateEquil_db(alphas_db(Y(2),aPR.WRT),betas_db(Y(2),aPR.WRT))-Y(6))/GateTimeCnst_db(alphas_db(Y(2),aPR.WRT),betas_db(Y(2),aPR.WRT));
dY(7) = (GateEquil_db(alphac_db(Y(2),aPR.WRT),betac_db(Y(2),aPR.WRT))-Y(7))/GateTimeCnst_db(alphac_db(Y(2),aPR.WRT),betac_db(Y(2),aPR.WRT));
dY(8) = (GateEquil_db(alphaq_db(Y(3)),betaq_db)-Y(8))/GateTimeCnst_db(alphaq_db(Y(3)),betaq_db);
end
function CaSatChi= Chi(Ca)
    CaSatChi = min(Ca/250,1);
end

function MInfsqr = MInfPR94(Vs,WRT)  
    MInfsqr = GateEquil_db(alpham_db(Vs,WRT),betam_db(Vs,WRT))^2;
end

end