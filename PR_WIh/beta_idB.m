function beta_i=beta_idB(v,Vihalf,wrt)
   v = v-(wrt+60);
   beta_i=exp(0.1581*(v-Vihalf));
end
