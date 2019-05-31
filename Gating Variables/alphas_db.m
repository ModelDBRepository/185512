function alpha_s =  alphas_db(Vd,WRT)    
    Vd = Vd-(WRT+60);
    num =1.6;
    den = 1 + exp(-0.072*(Vd-65));
    alpha_s = num./den;
end