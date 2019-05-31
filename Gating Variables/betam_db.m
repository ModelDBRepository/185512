function beta_m = betam_db(Vs,WRT)
%     n=size(Vs,1);
%     beta_m=zeros(n,1);
%     for i=1:n
%          tmpVs = Vs(i)-(WRT+60);
%          tmpnum = 0.28*(Vs(i)-40.1);
%          tmpden = exp((Vs(i)-40.1)/5)-1;
%          beta_m(i) = tmpnum/tmpden;
%     end;
     Vs=Vs-(WRT+60);
     num=0.28*(Vs-40.1);
     denexp=exp((Vs-40.1)/5);
     den=denexp-1;
     beta_m=num./den;
end