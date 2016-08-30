function [ T2, T1, T0 ] = cfdmat( gridx )
%CDFMAT:
% returns the matrices of second-order centered FD discretization
%         of u_xx, u_x and u

%INPUT:
%gridx: non_uniform grid on interval (a,b)

%OUTPUT:
%T2: discretization matrix of u_xx
%T1: discretization matrix of u_x
%T0: discretization matrix of u

    n = length(gridx) - 1; %number of grid
    m = n - 1; %number of interior points
    hx = (gridx(2:end) - gridx(1:end-1))'; %step size of non-uniform gridx
    
    hhd = hx(1:n-1).*hx(2:n);
    hhl = [hx(2:n-1).*(hx(2:n-1)+hx(3:n));1];
    hhu = [1;hx(2:n-1).*(hx(1:n-2)+hx(2:n-1))];
    
    T2 = spdiags([2./hhl -2./hhd 2./hhu], [-1 0 1], m, m);
    T1 = spdiags([-[hx(3:n);1]./hhl (hx(2:n)-hx(1:n-1))./hhd ...
        [1;hx(1:n-2)]./hhu], [-1 0 1], m, m);
    T0 = speye(m);
end