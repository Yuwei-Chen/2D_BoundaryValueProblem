% function [rhs, coefu, coefux, coefuy] = bc2(x, y, ib)
%
% returns the values of the BC coefficient functions and right side on (x, y)
% ib denotes the part of the boundary (x, y) belongs to
% ib = 1 on x = bx boundary line -- un = ux
% ib = 2 on y = ay boundary line -- un =-uy            7   4   8
% ib = 3 on x = ax boundary line -- un =-ux            3       1
% ib = 4 on y = by boundary line -- un = uy            6   2   5
% ib = 5 on (bx, ay) corner
% ib = 6 on (ax, ay) corner
% ib = 7 on (ax, by) corner
% ib = 8 on (bx, by) corner

function [rhs, coefu, coefux, coefuy] = bc2(x, y, ib)

coefu  =  1;
coefux =  0;
coefuy =  0;
% normal derivative?
if ib == 1
    coefux =  1;
    coefu  =  0;
%elseif ib == 2
%    coefuy = -1;
elseif ib == 3
    coefux = -1;
    coefu  =  0;
%elseif ib == 4
%    coefuy =  1;
elseif ib == 5
    coefux = 1;
    coefu  =  0;
%   coefuy = -1;
elseif ib == 6
    coefux = -1;
    coefu  =  0;
%   coefuy = -1;
elseif ib == 7
    coefux = -1;
    coefu  =  0;
%   coefuy = 1;
elseif ib == 8
    coefux = 1;
    coefu  =  0;
%   coefuy = 1;
end
coefu  =  1;
coefux =  0;
coefuy =  0;

%if x > 0.5, coefu = -2;, end; % discontinuity

[true1, true2, true3, true4, true5, true6] = truevd2(x, y);

rhs = coefu*true1 + coefux*true2 + coefuy*true4;
