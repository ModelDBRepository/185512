function alpha_c =  alphac_db(Vd,WRT)  

    %Issue with usinf if statements on arrays
    
    Vd = Vd-(WRT+60);
    alpha_c=zeros(1,size(Vd,2));
    for i=1:size(Vd,2)
    if Vd(i) <= 50       
        num = exp((Vd(i)-10)/11-(Vd(i)-6.5)/27);
        den = 18.975;
        alpha_c(i)=num/den;
    else
        alpha_c(i) = 2*exp((6.5-Vd(i))/27);
    end
    end
end