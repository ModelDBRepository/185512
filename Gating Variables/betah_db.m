function beta_h = betah_db(Vs,WRT)
    Vs = Vs-(WRT+60);
    num = 4;
    den = 1 + exp((40-Vs)/5);
    beta_h = num./den;
end