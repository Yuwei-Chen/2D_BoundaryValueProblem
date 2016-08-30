function [rhs, coefu, coefux, coefuxx, coefuy, coefuxy, coefuyy] = PDEcoefs(x, y)
%coefficients returns the values of the PDE coefficient functions
%             and RHS function g(x, y)
%   DEs o(x,y)*u_xx + p(x,y)*u_x + q(x,y)*u_yy + r(x,y)*u_y 
%   + s(x,y)*u_xy + t(x,y)*u = g(x,y), Dirichlet BCs

coefu   = -(x+y);
coefux  = sin(x)*sin(y);
coefuxx = exp(2*x+y);    % 1 + x*y; % 2^(0.5); % exp(x+y);
coefuy  = exp(-x-y);     % x*x + y*y;
coefuxy = x*y;           % x*y/5;
coefuyy = 1 + exp(x+y);  % 2 + x*y;
%%if x > 0.5, coefu = -2;, end;
%coefu   =  0; % -exp(y);
%coefux  = 1e1;
%coefuxx =  1; % 2+1/(1+y);
%coefuy  = 1e4; % 1+sin(y);
%coefuxy = .0;
%coefuyy =  1; % 3+y*y;

% rhs used to test when having the true value of functions
[true1, true2, true3, true4, true5, true6] = truevd2(x, y);
rhs = coefu*true1  + coefux*true2  + coefuxx*true3 ...
    + coefuy*true4 + coefuxy*true5 + coefuyy*true6;
%rhs = rhsfunc(x); %self-defined RHS function of PDE
end