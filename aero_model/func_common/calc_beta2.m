function beta = calc_beta2(u,v)
    %%              ^ u
    %   beta(-1,1)  |
    %    = 1/4pi    |     beta(1,1) = 7pi/4
    %               |
    %    ----------------------->v
    %    beta(-1,-1)|
    %    = 3/4pi    |     beta(1,-1) = 5pi/4
    %               |
    %%
    beta = abs(atan(v./u));   
    
    beta(intersect(find(u<0),find(v<0))) = pi-beta(intersect(find(u<0),find(v<0)));
    beta(intersect(find(u<0),find(v>0))) = pi+beta(intersect(find(u<0),find(v>0)));
    beta(intersect(find(u>0),find(v>0))) = 2*pi-beta(intersect(find(u>0),find(v>0)));
         
	beta((u&v) == 0) = 0;
    
end