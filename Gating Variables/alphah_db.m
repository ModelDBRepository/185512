function alpha_h =  alphah_db(Vs,WRT)    
    Vs = Vs-(WRT+60);
    alpha_h = 0.128*exp((17-Vs)/18);
end