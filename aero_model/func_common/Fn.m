function [A_F] = Fn(x,n,bias,U,alpha,beta)
% regressor generator 
% nth order Fourier series

if(nargin < 2)
    n=1;
end
if(nargin < 3)
    bias = 1;
end
if(nargin < 4)
    U = ones(size(x));    
end

if(nargin < 5)
    alpha = [];    
end

if(nargin < 6)
    beta = [];    
end

if bias == 0
    A_F = [];
else
    A_F = [ones(size(x)).*U];
end
A = cell(n,2);
if isempty(alpha)
    for i = 1:n
        A{i,1} = sin(i*x).*U;
        A{i,2} = cos(i*x).*U;
        A_F = [A_F A{i,1} A{i,2}];
    end
elseif isempty(beta)
    for i = 1:n
        A{i,1} = P1n(alpha,3,1,sin(i*x).*U);
        A{i,2} = P1n(alpha,3,1,cos(i*x).*U);
        A_F = [A_F A{i,1} A{i,2}];
    end
else
    for i = 1:n
        A{i,1} = P32(alpha,beta,sin(i*x).*U);
        A{i,2} = P32(alpha,beta,cos(i*x).*U);
        A_F = [A_F A{i,1} A{i,2}];
    end    
end
    
end

