function [A_p33] = P32(x1,x2,U)
% regressor generator 
% 3rd order polynomial function with 2 inputs
A1 = ones(size(x1)).*U;
A2 = x1.*U;
A3 = x2.*U;
A4 = x1.^2.*U;
A5 = x2.^2.*U;
A6 = x1.*x2.*U;
A7 = x1.^3.*U;
A8 = x2.^3.*U;
A9 = x1.*x2.^2.*U;
A10 = x1.^2.*x2.*U;

A_p33 = [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10];
end