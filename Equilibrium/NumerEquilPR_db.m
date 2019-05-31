function [numSSPR,diffProjFullEq,Jacob,eigJacob,nzeig] = NumerEquilPR_db(aPR,guessVsVd,VdsOut)

%% What does it do?
% This routine starts with a guess of the equillibrium Vs and Vd which
% might be found at VdsOut=0 mv or just really guessed. It then uses the
% equilibrium form of the PR equations which reduce to 2-d in Vs, and Vd
% It then uses MATLABs fsolve in 2d to get an equilibrium Vs and Vd and
% also calculates the corresponding values for the other 6 variables. It
% then goes on to compute the equilibrium a second way by using the 8
% equilibrium values found to initiate a search from the equilibrium of the
% full 8-d problem. The differences should be minimal and is output.

Vs= guessVsVd(1);
Vd= guessVsVd(2);

[ProjQuasiEqVsVd,quasiSSVsVdCahnscq,fval,exitflag,output,Jacob] = ProjGateEquilVsVd_db(aPR,Vs,Vd,VdsOut);

[numSSPR,diffProjFullEq,fval,exitflag,output,Jacob,eigJacob,nzeig] = GetFzeroSSPR_db(aPR,VdsOut,quasiSSVsVdCahnscq);

end