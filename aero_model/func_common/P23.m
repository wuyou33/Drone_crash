function [A_p23] = P23(x1,x2,x3)
% regressor generator 
% 2rd order polynomial function with 3 inputs

A1 = ones(size(x1));
A2 = x1;
A3 = x2;
A4 = x3;
A5 = x1.^2;
A6 = x2.^2;
A7 = x3.^3;
A8 = x1.*x2;
A9 = x2.*x3;
A10 = x3.*x1;

A_p23 = [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10];
end
