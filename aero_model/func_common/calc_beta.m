function beta = calc_beta(u,v)
    %%              ^ u
    %   beta(-1,1)  |
    %    = -1/4pi   |     beta(1,1) = pi/4
    %               |
    %    ----------------------->v
    %    beta(-1,-1)|
    %    = -3/4pi   |     beta(1,-1) = 3pi/4
    %               |
    %%
    beta = asin(v./sqrt(v.^2+u.^2));    
    beta(u<0) = sign(v(u<0)).*pi - beta(u<0);   
	beta((u&v) == 0) = 0;
    
end