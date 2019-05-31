 function alpha_m =  alpham_db(Vs,WRT)    
%     n=size(Vs,1);
%     alpha_m=zeros(n,1);
%     for i=1:n
%         tmpVs = Vs(i)-(WRT+60);
%         tmpnum = 0.32*(13.1-Vs(i));
%         tmpden = exp((13.1-Vs(i))/4)-1;
%         alpha_m(i) = tmpnum/tmpden;
%     end;
       Vs=Vs-(WRT+60);
     %  num=13.1*0.32-Vs*0.32;
       num=0.32*(13.1-Vs);
       denexp=(13.1-Vs)/4;
       den=exp(denexp)-1;
       alpha_m=num./den;
end