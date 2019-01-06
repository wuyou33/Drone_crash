function [A_p1n] = Fouriern(x1,n,bias,U)
% regressor generator 
% nrd order polynomial function with 1 inputs

if(nargin < 2)
    n=2;
end
if(nargin < 3)
    bias = 1;
end
if(nargin < 4)
    U = ones(size(x1));    
end


if bias == 0
    A_p1n = [];
else
    A_p1n = ones(size(x1));
end
A = cell(n,1);
for i = 1:n
    A{i} = x1.^i.*U;
    A_p1n = [A_p1n A{i}];
end

end
