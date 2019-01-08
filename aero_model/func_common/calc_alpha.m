function alpha = calc_alpha(u,w)
    %%              ^ w
    %   a(-1,1)     |
    %    = 3/4pi    |     alpha(1,1) = pi/4
    %               |
    %    ----------------------->u
    %    a(-1,-1)   |
    %    = -3/4pi   |     alpha(1,-1) = -pi/4
    %               |
    %%
    alpha = sign(w).*abs(atan(w./u));    
    alpha(u<0) = sign(w(u<0)).*pi - alpha(u<0);   
	alpha((u&w) == 0) = 0;
    
end