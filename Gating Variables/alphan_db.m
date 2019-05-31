function alpha_n =  alphan_db(Vs,WRT)    
    Vs = Vs-(WRT+60);
    num = 0.016*(35.1-Vs);
    den = exp((35.1-Vs)/5)-1;
    alpha_n = num./den;
end