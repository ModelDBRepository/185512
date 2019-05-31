function beta_c =  betac_db(Vd,WRT)    
    Vd = Vd-(WRT+60);
    beta_c=zeros(1,size(Vd,2));
    for i=1:size(Vd,2)
    if Vd(i) <= 50       
        beta_c(i) = 2*exp((6.5-Vd(i))/27) - alphac_db(Vd(i),WRT);
    else
        beta_c(i) = 0;
    end
    end
end