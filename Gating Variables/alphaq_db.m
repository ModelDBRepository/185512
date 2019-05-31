function alpha_q =  alphaq_db(Ca)    
    alpha_q = min(0.00002*Ca,0.01);
end