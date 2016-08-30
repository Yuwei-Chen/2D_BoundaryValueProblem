function [ A ] = cfdmat2d( gridx, gridy, coefs )
%CFDMAT2D:
% returns the matrix of second-order centered FD discretization
%         of % of 2-dimensional second-order DEs o(x,y)*u_xx + p(x,y)*u_x 
%         + q(x,y)*u_yy + r(x,y)*u_y + s(x,y)*u_xy + t(x,y)*u = g(x,y)
%         Dirichlet BCs

%INPUT:
%gridx: non_uniform grid on interval (ax,bx)
%gridy: non_uniform grid on interval (ay,by)
%coefs: values of o, p, q, r, s, t

%OUTPUT:
%A: discretization matrix of 2-dimensional DEs

[T2x, T1x, T0x] = cfdmat(gridx);
[T2y, T1y, T0y] = cfdmat(gridy);

A = spdiag(coefs(:,1))*kron(T0x, T0y)...
  + spdiag(coefs(:,2))*kron(T1x, T0y)...
  + spdiag(coefs(:,3))*kron(T2x, T0y)...
  + spdiag(coefs(:,4))*kron(T0x, T1y)...
  + spdiag(coefs(:,5))*kron(T1x, T1y)...
  + spdiag(coefs(:,6))*kron(T0x, T2y);
end