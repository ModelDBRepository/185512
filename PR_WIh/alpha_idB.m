function alpha_i=alpha_idB(v,Vihalf,wrt)
   v = v-(wrt+60);
   alpha_i=exp(-0.1054*(v-Vihalf));
end