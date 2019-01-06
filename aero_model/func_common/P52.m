function [A_p52] = P52(x1,x2,U)
% regressor generator 
% 5rd order polynomial function with 2 inputs
if(nargin < 3)
    U = ones(size(x1));
end
A1 = ones(size(x1)).*U;
A2 = x1.*U;
A3 = x2.*U;
A4 = x2.^2.*U;
A5 = x1.*x2.*U;
A6 = x2.^2.*U;
A7 = x1.^3.*U;
A8 = x1.^2.*x2.*U;
A9 = x1.*x2.^2*U;
A10 = x2.^3.*U;
A11 = x1.^4.*U;
A12 = x1.^3.*x2.*U;
A13 = x1.^2.*x2.^2.*U;
A14 = x1.*x2.^3.*U;
A15 = x2.^4.*U;
A16 = x1.^5.*U;
A17 = x1.^4.*x2.*U;
A18 = x1.^3.*x2.^2.*U;
A19 = x1.^2.*x2.^3.*U;
A20 = x1.*x2.^4.*U;
A21 = x2.^5.*U;

A_p52 = [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,...
         A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21];
end