function [CT, T] = generate_CT_T(gmm)
while true
    while true
        R = round(random(gmm));


        if R(1) >= 0 && R(1)<=6 && R(2)>=0 && R(2)<=40
           break; 
        end
    end
    
    if (R(1)==0||R(1)==1) && R(2)==0
        break;
    end
    
    if (R(1)==0||R(1)==1) && R(2) > 0
        continue
    end
    
break;
end
    
CT = R(1);

T = R(2);
end