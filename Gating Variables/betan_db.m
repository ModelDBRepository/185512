function beta_n = betan_db(Vs,WRT)
    Vs = Vs-(WRT+60);
    beta_n = 0.25*exp(0.5-0.025*Vs);
end