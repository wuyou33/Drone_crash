function B = calc_Bs(X,U,P)

v = X(2); 
r = X(3);
% a = X(4);
% q = X(5);
% w = X(6);
% u = U(1);
s = P(4);

B =[0;0;0;0;0;
 -s/(v*sin(r))];

end