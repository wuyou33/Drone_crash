function [A_p33] = P33(x1,x2,x3,U)
% regressor generator 
% 3rd order polynomial function with 3 inputs
if (nargin < 3)
    x3 = ones(size(x1));
end
if(nargin < 4)
    U = ones(size(x1));
end
A1 = ones(size(x1)).*U;
A2 = x1.*U;
A3 = x2.*U;
A4 = x3.*U;
A5 = x1.^2.*U;
A6 = x2.^2.*U;
A7 = x3.^2.*U;
A8 = x1.*x2.*U;
A9 = x2.*x3.*U;
A10 = x3.*x1.*U;
A11 = x1.^3.*U;
A12 = x2.^3.*U;
A13 = x3.^3.*U;
A14 = x1.*x2.^2.*U;
A15 = x1.*x3.^2.*U;
A16 = x1.^2.*x2.*U;
A17 = x1.^2.*x3.*U;
A18 = x2.^2.*x3.*U;
A19 = x2.*x3.^2.*U;
A20 = x1.*x2.*x3.*U;

A_p33 = [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,...
         A11,A12,A13,A14,A15,A16,A17,A18,A19,A20];
end