function beta_s = betas_db(Vd,WRT)
    Vd = Vd-(WRT+60);
    num = 0.02*(Vd-51.1);
    den = exp((Vd-51.1)/5)-1;
    beta_s = num./den;
end