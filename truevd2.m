% function [true1, true2, true3, true4, true5, true6] = U(x, y)
%
% returns the values of the solution function,
% and derivatives up to 2nd order on x, y, i.e. u, ux, uxx, uy, uxy, uyy
% x can be a vector/matrix. If x is a vector/matrix,
% then y must be either a scalar or a vector/matrix of same size as x.

function [true1, true2, true3, true4, true5, true6] = truevd2(x, y)

global Uno Uname;
global RR SS QQ KK AA BB EE TT ax bx;
global D4X D3X; % temporary and incomplete!!!
global W1 W2 ax bx gam eta kap ee kap3 kap2 kap1 kap0;

oo = ones(size(x)); zz = zeros(size(x));
yy = y;
if length(y) == 1, y = repmat(y, size(x));, end
switch Uno
case {700, 702}
    Uname = 'European Call (not all derivatives computed)';
    S = x; t = y; r = RR; s = SS; q = QQ; E = EE; T = TT;
    if (t == T)
        true1 = max(S-E, zz); %if (S >= E*oo), true2 = oo; else true2 = zz; end
        for i=1:length(x)
            if (S(i) >= EE), true2(i) = 1; else true2(i) = 0; end
        end
        true3 = zz; true4 = zz; true5 = zz; true6 = zz;
        return
    end
    zi = find(x <= 0); x(zi) = 10^(-20); S = x;
    d1 = (log(S/E) + (r - q + s^2/2)*(T-t))./(s*sqrt(T-t));
    d2 = (log(S/E) + (r - q - s^2/2)*(T-t))./(s*sqrt(T-t));
    Nd1 = mynormcdf(d1);
    Nd2 = mynormcdf(d2);
    true1 = S.*Nd1 - E*exp(-r*(T-t)).*Nd2;
    %true2 = Nd1;
    true1(zi) = zeros(size(zi)); %true2(zi) = ones(size(zi));
    %true3 = zz; true4 = zz; true5 = zz; true6 = zz;
    [true2, true3, true4] = egreeks(S, E, r, s, q, T, T-t, 'c');

case {701, 703}
    Uname = 'European Put (not all derivatives computed)';
%erf(x) = 2/sqrt(pi) * integral from 0 to x of exp(-t^2) dt
    S = x; t = y; r = RR; s = SS; q = QQ; E = EE; T = TT;
    if (t == T)
        true1 = max(E-S, zz); %if (S <= E*oo), true2 = -oo; else true2 = zz; end
        for i=1:length(x)
            if (S(i) <= EE), true2(i) = -1; else true2(i) = 0; end
        end
        true3 = zz; true4 = zz; true5 = zz; true6 = zz;
	return
    end
    zi = find(x <= 0); x(zi) = 1e-20; S = x;
    d1 = (log(S/E) + (r - q + s^2/2)*(T-t))./(s*sqrt(T-t));
    d2 = (log(S/E) + (r - q - s^2/2)*(T-t))./(s*sqrt(T-t));
    Nd1 = mynormcdf(-d1);
    Nd2 = mynormcdf(-d2);
    true1 = E*exp(-(r-q)*(T-t)).*Nd2 - S.*Nd1;
    %true2 = -Nd1;
    %true1(zi) = E*exp(-r*(T-t))*ones(size(zi)); %true2(zi) = zeros(size(zi));
    %true3 = zz; true4 = zz; true5 = zz; true6 = zz;
    [true2, true3, true4] = egreeks(S, E, r, s, q, T, T-t, 'p');

case {704}
    Uname = 'European Call forward (not all derivatives computed)';
    S = x; t = y; r = RR; s = SS; q = QQ; E = EE; T = TT;
    if (t == 0)
        true1 = max(S-E, zz); %if (S >= E*oo), true2 = oo; else true2 = zz; end
        for i=1:length(x)
            if (S(i) >= EE), true2(i) = 1; else true2(i) = 0; end
        end
        true3 = zz; true4 = zz; true5 = zz; true6 = zz;
        return
    end
    zi = find(x <= 0); x(zi) = 10^(-20); S = x;
    d1 = (log(S/E) + (r - q + s^2/2)*(t))./(s*sqrt(t));
    d2 = (log(S/E) + (r - q - s^2/2)*(t))./(s*sqrt(t));
    Nd1 = normcdf(d1);
    Nd2 = normcdf(d2);
    true1 = S.*Nd1 - E*exp(-r*(t)).*Nd2;
    %true2 = Nd1;
    true1(zi) = zeros(size(zi)); %true2(zi) = ones(size(zi));
    %true3 = zz; true4 = zz; true5 = zz; true6 = zz;
    [true2, true3, true4] = egreeks(S, E, r, s, q, T, t, 'c');

case {705}
    Uname = 'European Put forward (not all derivatives computed)';
%erf(x) = 2/sqrt(pi) * integral from 0 to x of exp(-t^2) dt
    S = x; t = y; r = RR; s = SS; q = QQ; E = EE; T = TT;
    if (t == 0)
        true1 = max(E-S, zz); %if (S <= E*oo), true2 = -oo; else true2 = zz; end
        for i=1:length(x)
            if (S(i) <= EE), true2(i) = -1; else true2(i) = 0; end
        end
        true3 = zz; true4 = zz; true5 = zz; true6 = zz;
	return
    end
    zi = find(x <= 0); x(zi) = 1e-20; S = x;
    d1 = (log(S/E) + (r - q + s^2/2)*(t))./(s*sqrt(t));
    d2 = (log(S/E) + (r - q - s^2/2)*(t))./(s*sqrt(t));
    Nd1 = mynormcdf(-d1);
    Nd2 = mynormcdf(-d2);
    true1 = E*exp(-(r-q)*(t)).*Nd2 - S.*Nd1;
    %true2 = -Nd1;
    true1(zi) = E*exp(-r*(t(1)))*ones(size(zi)); %true2(zi) = zeros(size(zi));
    %true3 = zz; true4 = zz; true5 = zz; true6 = zz;
    [true2, true3, true4] = egreeks(S, E, r, s, q, T, t, 'p');

case {7010}
    Uname = 'American Put forward (computed at start and at boundary)';
    S = x; t = y; r = RR; q = QQ; s = SS; E = EE; T = TT;
    if (t == 0)
        true1 = max(E-S, zz);
        for i=1:length(x)
            if (S(i) < E), true2(i) = -1; else true2(i) = 0; end
        end
        true3 = zz; true4 = zz; true5 = zz; true6 = zz;
        return
    end
    if (S == bx)
        true1 = zz; true2 = zz; true3 = zz; true4 = zz; true5 = zz; true6 = zz;
        return
    end
    if (S == ax)
        true1 = E;  true2 = -1; true3 = zz; true4 = zz; true5 = zz; true6 = zz;
        return
    end

case {706}
    Uname = 'Jump diffusion, no derivative computed';
    r = RR; T = TT; K = EE; sigma = SS; t =yy;
    S = x;
    global lambda type tolMerton jumpDensity;
    % type = 1 for call, type = 2 for put
    if (abs(T-t) < 1e-9)
        if (type == 1)
            true1 = max(S-K, 0);
            if (S > K*oo), true2 = oo; else true2 = zz; end
        elseif (type == 2)
            true1 = max(K-S, 0);
            if (S < K*oo), true2 = -oo; else true2 = zz; end
        end
        true3 = zz; true4 = zz; true5 = zz; true6 = zz;
        return;
    end
    true1 = zeros(size(x));
    if (jumpDensity == 1)
        global gamma mu
        for i = 1:length(x)
            true1(i) = merton(x(i), K, r, T-t, sigma, gamma, mu, lambda, ...
                              type, tolMerton);
        end
    elseif (jumpDensity == 2)
        global p eta1 eta2
        true1 = doubexp(x, K, r, T-t, sigma, lambda, p, eta1, eta2, ...
                        type, tolMerton);
    end
    true2 = zz; true3 = zz; true4 = zz; true5 = zz; true6 = zz;

case {7061}
    Uname = 'Jump diffusion, up and out call';
    r = RR; T = TT; K = EE; sigma = SS; t =yy;
    S = x;
    global lambda type tolMerton jumpDensity barrier bar_val;
    % type = 1 for call, type = 2 for put
    if (abs(T-t) < 1e-9)
        true1 = zeros(length(S), 1);
        if (type == 1)
            for i = 1:length(S)
                true1(i) = max(S(i)-K, 0);
                if (S(i) > K*oo), true2(i) = 1; else true2(i) = 0; end
                if (S(i) >= bar_val) true1(i) = 0; end
            end
        end
        true3 = zz; true4 = zz; true5 = zz; true6 = zz;
        return;
    end

    true1 = zeros(size(x));

    if (lambda)
        % incorrect
        if (jumpDensity == 1)
            global gamma mu
            for i = 1:length(x)
                true1(i) = merton(x(i), K, r, T-t, sigma, gamma, mu, lambda, ...
                                  type, tolMerton);
            end
        elseif (jumpDensity == 2)
            global p eta1 eta2
            true1 = doubexp(x, K, r, T-t, sigma, lambda, p, eta1, eta2, ...
                            type, tolMerton);
        end
    else
        B = bar_val;
        for i = 1:length(x)
            if (x(i) == 0) true1(i) = 0; continue; end
            if (x(i) >= B) true1(i) = 0 ; continue; end
            true1(i) = barrier_price(x(i), K, B, T-t, sigma, r);
        end
    end

    true2 = zz; true3 = zz; true4 = zz; true5 = zz; true6 = zz;

case {185}
    Uname = 'exp(-a*((x-c)^2+(y-d)^2))*(x^2-x)*(y^2-y)';
    a = 1000; c = .8; d = .8;
    ee = exp(-a*((x-c).^2+(y-d).^2));
    yy = (y.^2-y);
    xx = (x.^2-x);
    xy = xx.*yy;
    x2 = (2.*x-2*c);
    y2 = (2.*y-2*d);
    true1 = ee.*xy;
    true2 = -a.*x2.*true1 + ee.*(2.*x-1).*yy;
    true3 = -2.*a.*true1 + ee.*(a^2.*x2.^2.*xy ...
            -2.*a.*x2.*(2.*x-1).*yy + 2.*yy);
    true4 = -a.*y2.*true1 + ee.*(2.*y-1).*xx;
    true5 =  a^2.*y2.*x2.*true1 + ee.*(-a.*x2.*xx.*(2.*y-1) ...
            -a.*y2.*(2.*x-1).*yy + (2.*x-1).*(2.*y-1));
    true6 = -2*a.*true1 + ee.*(a^2.*y2.^2.*xy ...
            -2.*a.*y2.*(2.*y-1).*xx + 2.*xx);
case {165}
    Uname = 'x.^(13/2) .* y.^(13/2)';
    true1 =          x.^(13/2) .*          y.^(13/2);
    true2 =  13/2 .* x.^(11/2) .*          y.^(13/2);
    true3 = 143/4 .* x.^(9/2)  .*          y.^(13/2);
    true4 =          x.^(13/2) .*  13/2 .* y.^(11/2);
    true5 =  13/2 .* x.^(11/2) .*  13/2 .* y.^(11/2);
    true6 =          x.^(13/2) .* 143/4 .* y.^(9/2);
case {163}
    Uname = 'x.^(13/2) .* exp(a*y) a = -1'; a = -1; ya = a*y;
    true1 =         x.^(13/2) .* exp(ya);
    true2 =  13/2 .*x.^(11/2) .* exp(ya);
    true3 = 143/4 .*x.^(9/2)  .* exp(ya);
    true4 =         x.^(13/2) .* exp(ya).*a;
    true5 =  13/2 .*x.^(11/2) .* exp(ya).*a;
    true6 =         x.^(13/2) .* exp(ya).*a.^2;
case {155}
    Uname = 'x.^(11/2) .* y.^(11/2)';
    true1 =         x.^(11/2) .*        y.^(11/2);
    true2 =  11/2 .*x.^(9/2)  .*        y.^(11/2);
    true3 =  99/4 .*x.^(7/2)  .*        y.^(11/2);
    true4 =         x.^(11/2) .* 11/2 .*y.^(9/2);
    true5 =  11/2 .*x.^(9/2)  .* 11/2 .*y.^(9/2);
    true6 =         x.^(11/2) .* 99/4 .*y.^(7/2);
case {153}
    Uname = 'x.^(11/2) .* exp(a*y) a = -1'; a = -1; ya = a*y;
    true1 =         x.^(11/2) .* exp(ya);
    true2 =  11/2 .*x.^(9/2)  .* exp(ya);
    true3 =  99/4 .*x.^(7/2)  .* exp(ya);
    true4 =         x.^(11/2) .* exp(ya).*a;
    true5 =  11/2 .*x.^(9/2)  .* exp(ya).*a;
    true6 =         x.^(11/2) .* exp(ya).*a.^2;
case {145}
    Uname ='x.^(9/2) .* y.^(9/2)';
    true1 =        x.^(9/2) .*        y.^(9/2);
    true2 =  9/2 .*x.^(7/2) .*        y.^(9/2);
    true3 = 63/4 .*x.^(5/2) .*        y.^(9/2);
    true4 =        x.^(9/2) .*  9/2 .*y.^(7/2);
    true5 =  9/2 .*x.^(7/2) .*  9/2 .*y.^(7/2);
    true6 =        x.^(9/2) .* 63/4 .*y.^(5/2);
case {143}
    Uname ='x.^(9/2) .* exp(a*y) a = -1'; a = -1; ya = a*y;
    true1 =        x.^(9/2) .* exp(ya);
    true2 =  9/2 .*x.^(7/2) .* exp(ya);
    true3 = 63/4 .*x.^(5/2) .* exp(ya);
    true4 =        x.^(9/2) .* exp(ya).*a;
    true5 =  9/2 .*x.^(7/2) .* exp(ya).*a;
    true6 =        x.^(9/2) .* exp(ya).*a.^2;
case {144}
    Uname = '(x.^(13/2) -2*x.^(11/2) +x.^(9/2)).*(y.^(13/2) -2*y.^(11/2) +y.^(9/2))';
    true1 =  (         x.^(13/2) -    2 .* x.^(11/2) +         x.^(9/2)) ...
          .* (         y.^(13/2) -    2 .* y.^(11/2) +         y.^(9/2));
    true2 =  ( 13/2 .* x.^(11/2) -   11 .* x.^(9/2)  +  9/2 .* x.^(7/2)) ...
          .* (         y.^(13/2) -    2 .* y.^(11/2) +         y.^(9/2));
    true3 =  (143/4 .* x.^(9/2)  - 99/2 .* x.^(7/2)  + 63/4 .* x.^(5/2)) ...
          .* (         y.^(13/2) -    2 .* y.^(11/2) +         y.^(9/2));
    true4 =  (         x.^(13/2) -    2 .* x.^(11/2) +         x.^(9/2)) ...
          .* ( 13/2 .* y.^(11/2) -   11 .* y.^(9/2)  +  9/2 .* y.^(7/2));
    true5 =  ( 13/2 .* x.^(11/2) -   11 .* x.^(9/2)  +  9/2 .* x.^(7/2)) ...
          .* ( 13/2 .* y.^(11/2) -   11 .* y.^(9/2)  +  9/2 .* y.^(7/2));
    true6 =  (         x.^(13/2) -    2 .* x.^(11/2) +         x.^(9/2)) ...
          .* (143/4 .* y.^(9/2)  - 99/2 .* y.^(7/2)  + 63/4 .* y.^(5/2));
case {135}
    Uname ='x.^(7/2) .* y.^(7/2)';
    true1 =        x.^(7/2) .*        y.^(7/2);
    true2 =  7/2 .*x.^(5/2) .*        y.^(7/2);
    true3 = 35/4 .*x.^(3/2) .*        y.^(7/2);
    true4 =        x.^(7/2) .*  7/2 .*y.^(5/2);
    true5 =  7/2 .*x.^(5/2) .*  7/2 .*y.^(5/2);
    true6 =        x.^(7/2) .* 35/4 .*y.^(3/2);
case {120}
    Uname ='exp(x)+exp(y)';
    true1 = exp(x)+exp(y);   true2 = exp(x);          true3 = exp(x);
    true4 = exp(y);          true5 = zz;              true6 = exp(y);
    D4X = exp(x);
case {119}
    Uname ='exp(x+y)';
    true1 = exp(x+y);        true2 = exp(x+y);        true3 = exp(x+y);
    true4 = exp(x+y);        true5 = exp(x+y);        true6 = exp(x+y);
    D4X = exp(x+y);
case {118}
    Uname ='(x.*x - x).*(y.*y - y).*exp(x+y)';
    true1 = (x.*x - x).*(y.*y - y).*exp(x+y);
    true2 = (2.*x - 1).*(y.*y - y).*exp(x+y) + true1;
    true3 = (2.*(y.*y - y) + (2.*x - 1).*(y.*y - y)).*exp(x+y) + true2;
    true4 = (x.*x - x).*(2.*y - 1).*exp(x+y) + true1;
    true5 = (2.*x - 1).*(2.*y - 1).*exp(x+y) ...
          + (2.*x - 1).*(y.*y - y).*exp(x+y) + true4;
    true6 = (2.*(x.*x - x) + (2.*y - 1).*(x.*x - x)).*exp(x+y) + true4;
case {117}
    Uname ='(x.*x - x).*(y.*y - y).*exp(-a*(x-x0).^2 -b*(y-y0).^2)';
    a = 10; b = 40; x0 = 0.7; y0 = 0.4;
    true1 = (x.*x - x).*(y.*y - y).*exp(-a*(x-x0).^2 -b*(y-y0).^2);
    u = exp(-a*(x-x0).^2 -b*(y-y0).^2);
    true2 = ((2*x-1) - 2*a*(x-x0).*(x.*x-x)).*(y.*y-y).*u;
    true3 = (2 -4*(2*x-1)*a*(x-x0) -2*(x.*x-x)*a + 4*(x.*x-x)*a^2.*(x-x0).^2).*(y.*y-y).*u;
    true4 = ((2*y-1) - 2*b*(y-y0).*(y.*y-y)).*(x.*x-x).*u;
    true5 = ((2*x-1).*(2*y-1) - 2*b*(2*x-1)*(y-y0).*(y.*y-y)).*u ...
          + (-2*a*(2*y-1)*(x-x0).*(x.*x-x) + 4*(x.*x - x).*(y.*y - y)*a*b*(x-x0).*(y-y0)).*u;
    true6 = (2 -4*(2*y-1)*b*(y-y0) -2*(y.*y-y)*b + 4*(y.*y-y)*b^2.*(y-y0).^2).*(x.*x-x).*u;
case {116}
    Uname ='(x.*x - x).*exp(x+y)';
    true1 = (x.*x - x).*exp(x+y);
    true2 = (2.*x - 1).*exp(x+y) + true1;
    true3 = (2 + (2.*x - 1)).*exp(x+y) + true2;
    true4 = true1;
    true5 = true2;
    true6 = true1;
    D4X = ((x.*x - x) + 4*(2*x-1) + 12).*exp(x+y);
case {115}
% u = 0, uxxxx = 0 on x = 0, 1
    Uname ='(x + x.^2 + 2*x.^3 - 6*x.^5 + 2*x.^6).*exp(y) -- u, uxxxx = 0 on x = 0, 1';
    true1 = (x + x.^2 + 2*x.^3 -  6*x.^5 + 2*x.^6).*exp(y);
    true2 = (1 +2*x   + 6*x.^2 - 30*x.^4 +12*x.^5).*exp(y);
    true3 = (    2    +12*x    -120*x.^3 +60*x.^4).*exp(y);
    true4 = true1;
    true5 = true2;
    true6 = true1;
    D4X   = (                  -720*x  + 720*x.^2).*exp(y);
case {114}
    Uname ='(x^4 - x).*exp(y)';
% u = 0 on x = 0, 1
    tx1 = x.^4 - x; tx2 = 4*x.^3 - 1; tx3 = 12*x.*x; ty = exp(y);
    true1 = tx1.*ty;         true2 = tx2.*ty;         true3 = tx3.*ty;
    true4 = tx1.*ty;         true5 = tx2.*ty;         true6 = tx1.*ty;
    D4X = 24*oo.*exp(y);
case {113}
    Uname ='sin(pi*x).*exp(y*a) a = -1'; px = pi*x; a = -1; ya = a*y;
    true1 =      sin(px).*exp(ya);
    true2 =   pi*cos(px).*exp(ya);
    true3 =-pi^2*sin(px).*exp(ya);
    true4 =      sin(px).*exp(ya).*a;
    true5 =   pi*cos(px).*exp(ya).*a;
    true6=       sin(px).*exp(ya).*a.^2;
case {112}
    Uname ='sin(2*pi*x).*exp(y)'; tp = 2*pi; tpx = tp*x;
    true1 = sin(tpx).*exp(y); true2 = tp*cos(tpx).*exp(y); true3=-tp^2*sin(tpx).*exp(y);
    true4 = sin(tpx).*exp(y); true5 = tp*cos(tpx).*exp(y); true6= sin(tpx).*exp(y);
case {111}
    Uname ='sin(x).*exp(y)';
    true1 = sin(x).*exp(y); true2 = cos(x).*exp(y); true3=-sin(x).*exp(y);
    true4 = sin(x).*exp(y); true5 = cos(x).*exp(y); true6= sin(x).*exp(y);
    D4X = true1;
case {110}
    Uname ='sin(x).*exp(-y)';
    true1 = sin(x).*exp(-y); true2 = cos(x).*exp(-y); true3=-sin(x).*exp(-y);
    true4 =-sin(x).*exp(-y); true5 =-cos(x).*exp(-y); true6= sin(x).*exp(-y);
    D4X = true1;
case {108}
    Uname ='sin(x).*sin(y)';
    true1 = sin(x).*sin(y);  true2 = cos(x).*sin(y);  true3 =-sin(x).*sin(y);
    true4 = sin(x).*cos(y);  true5 = cos(x).*cos(y);  true6 =-sin(x).*sin(y);
    D4X = true1;
case {107}
    Uname ='sin(x+y) + 4';
    true1 = sin(x+y) + 4;    true2 = cos(x+y);        true3 =-sin(x+y);
    true4 = cos(x+y);        true5 =-sin(x+y);        true6 =-sin(x+y);
    D4X = true1;
case {106}
    Uname ='cos(x+y-pi)';
    true1 = cos(x+y-pi);     true2 =-sin(x+y-pi);     true3 =-cos(x+y-pi);
    true4 =-sin(x+y-pi);     true5 =-cos(x+y-pi);     true6 =-cos(x+y-pi);
    D4X = true1;
case {105}
    Uname ='cos(x).*y';
    true1 = cos(x).*y;       true2 =-sin(x).*y;       true3 =-cos(x).*y;
    true4 = cos(x);          true5 =-sin(x);          true6 = zz;
    D4X = true1;

case {104, -104}
    Uname ='exp(-y)*kap*exp(-x.^2/eta^2) + terms for discontinuity of derv u at W1, W2';
%kap = exp(-y);
true1 = kap*exp(-x.^2/eta^2) .* exp(-y);
true2 =-2*x/eta^2*kap.*exp(-x.^2/eta^2) .* exp(-y);
true3 =(-2/eta^2*kap*exp(-x.^2/eta^2) + 4*x.^2/eta^4*kap.*exp(-x.^2/eta^2)) .* exp(-y);
n = length(x);
for i = 1:n
    if (ax <= x(i)) & (x(i) <= W1)
    elseif (W1 < x(i)) & (x(i) < W2)
        true1(i) = true1(i) + (x(i)-W1).*(x(i)-W2) .* exp(-y(i));
        true2(i) = true2(i) + (2*x(i) - W1 - W2) .* exp(-y(i));
        true3(i) = true3(i) + 2*exp(-y(i));
    elseif (W2 <= x(i)) & (x(i) <= bx)
    end
end
true4 = -true1;
true5 = -true2;
true6 = true1;
D3X = (12*x/eta^4 - 8*x.^3/eta^6).*true1;
if (size(D3X, 1) == 1), D3X = D3X';, end
D4X = (12/eta^4 - 48*x.^2/eta^6 + 16*x.^4/eta^8).*true1;
if (size(D4X, 1) == 1), D4X = D4X';, end

case {104.1}
    Uname ='kap*exp(-x.^2/eta^2) + terms for discontinuity of derv u at W1, W2';
true1 = kap*exp(-x.^2/eta^2);
true2 =-2*x/eta^2*kap.*exp(-x.^2/eta^2);
true3 =-2/eta^2*kap.*exp(-x.^2/eta^2) + 4*x.^2/eta^4*kap.*exp(-x.^2/eta^2);
n = length(x);
for i = 1:n
    if (ax <= x(i)) & (x(i) <= W1)
    elseif (W1 < x(i)) & (x(i) < W2)
        true1(i) = true1(i) + (x(i)-W1).*(x(i)-W2);
        true2(i) = true2(i) + (2*x(i) - W1 - W2);
        true3(i) = true3(i) + 2;
    elseif (W2 <= x(i)) & (x(i) <= bx)
    end
end
true4 = zz;
true5 = zz;
true6 = zz;
D3X = (12*x/eta^4 - 8*x.^3/eta^6).*true1;
if (size(D3X, 1) == 1), D3X = D3X';, end
D4X = (12/eta^4 - 48*x.^2/eta^6 + 16*x.^4/eta^8).*true1;
if (size(D4X, 1) == 1), D4X = D4X';, end

case {103.5, -103.5}
    Uname ='exp(-y)*(kap3*x.^3 + kap*x.^2) + terms for discontinuity of derv u at W1, W2';
true1 = (kap3*x.^3 + kap*x.^2 + kap1*x).*exp(-y);
true2 = (3*kap3*x.^2 + 2*kap*x + kap1).*exp(-y);
true3 = (6*kap3*x + 2*kap).*exp(-y);
n = length(x);
for i = 1:n
if (ax <= x(i)) & (x(i) <= W1)
elseif (W1 < x(i)) & (x(i) < W2)
    true1(i) = true1(i) + (x(i)-W1).*(x(i)-W2) .* exp(-y(i));
    true2(i) = true2(i) + (2*x(i) - W1 - W2) .* exp(-y(i));
    true3(i) = true3(i) + 2 .* exp(-y(i));
elseif (W2 <= x(i)) & (x(i) <= bx)
end
end
true4 =-true1;
true5 =-true2;
true6 = true1;

case {103.3}
    Uname ='kap3*x.^3 + kap*x.^2 + terms for discontinuity of derv u at W1, W2';
true1 = kap3*x.^3 + kap*x.^2 + kap1*x;
true2 = 3*kap3*x.^2 + 2*kap*x + kap1;
true3 = 6*kap3*x + 2*kap;
n = length(x);
for i = 1:n
if (ax <= x(i)) & (x(i) <= W1)
elseif (W1 < x(i)) & (x(i) < W2)
    true1(i) = true1(i) + (x(i)-W1).*(x(i)-W2);
    true2(i) = true2(i) + 2*x(i) - W1 - W2;
    true3(i) = true3(i) + 2;
elseif (W2 <= x(i)) & (x(i) <= bx)
end
end
true4 = zz;
true5 = zz;
true6 = zz;

case {103}
    Uname ='exp(-y)*(kap*x.^2 + terms for discontinuity of derv u at W1, W2)';
true1 = kap*x.^2 .* exp(-y);
true2 = 2*kap*x .* exp(-y);
true3 = 2*kap*oo .* exp(-y);
n = length(x);
for i = 1:n
    if (ax <= x(i)) & (x(i) <= W1)
    elseif (W1 < x(i)) & (x(i) < W2)
        true1(i) = true1(i) + (x(i)-W1).*(x(i)-W2) .* exp(-y(i));
        true2(i) = true2(i) + (2*x(i) - W1 - W2) .* exp(-y(i));
        true3(i) = true3(i) + 2 .* exp(-y(i));
    elseif (W2 <= x(i)) & (x(i) <= bx)
    end
end
true4 = -true1;
true5 = -true2;
true6 =  true1;

case {103.1}
    Uname ='kap*x.^2 + terms for discontinuity of derv u at W1, W2';
true1 = kap*x.^2;
true2 = 2*kap*x;
true3 = 2*kap*oo;
n = length(x);
for i = 1:n
    if (ax <= x(i)) & (x(i) <= W1)
    elseif (W1 < x(i)) & (x(i) < W2)
        true1(i) = true1(i) + (x(i)-W1).*(x(i)-W2);
        true2(i) = true2(i) + (2*x(i) - W1 - W2);
        true3(i) = true3(i) + 2;
    elseif (W2 <= x(i)) & (x(i) <= bx)
    end
end 
true4 = zz;
true5 = zz;
true6 = zz;

case {-103.1}
    Uname ='initial & boundary data only - piece. lin. schewed (hat) with DD at W1, W2';
% point of peak (1*W1+2*W2)/3, value at peak = 1; slopes 3/4, -3/2
ppeak = (1*W1+2*W2)/3; vpeak = 1;
slp1 = vpeak/(ppeak-W1); slp2 = vpeak/(ppeak-W2);
true1 = zz;
true2 = zz;
true3 = zz;
n = length(x);
for i = 1:n
    if (ax <= x(i)) & (x(i) <= W1)
    elseif (W1 < x(i)) & (x(i) <= ppeak)
        true1(i) = slp1*x(i);
        true2(i) = slp1;
        true3(i) = 0;
    elseif (ppeak < x(i)) & (x(i) < W2)
        true1(i) = vpeak+slp2*(x(i)-ppeak);
        true2(i) = slp2;
        true3(i) = 0;
    elseif (W2 <= x(i)) & (x(i) <= bx)
    end
end
true4 = zz;
true5 = zz;
true6 = zz;

case {-103.2}
    Uname ='initial & boundary data only - piece. quadr. with DD at -1, +1';
true1 = zz;
true2 = zz;
true3 = zz;
n = length(x);
for i = 1:n
    if (ax <= x(i)) & (x(i) <= -1)
    elseif (-1 < x(i)) & (x(i) < 1)
        true1(i) = -x(i)^2 + 1;
        true2(i) = -2*x(i);
        true3(i) = -2;
    elseif (1 <= x(i)) & (x(i) <= bx)
    end
end
true4 = zz;
true5 = zz;
true6 = zz;

case {-103.3}
    Uname ='initial & boundary data only - mix of exp and cub with DD at 0, 2';
% matching up to including second derivative at 7/4
% discontinuous first derivative at 2, fully conditnuous elsewhere
true1 = kap*exp(-(x-1.75).^2/eta^2);
true2 = kap*exp(-(x-1.75).^2/eta^2) .* (-2*(x-1.75))/eta^2;
true3 = kap*exp(-(x-1.75).^2/eta^2) .* ( 4*(x-1.75).^2)/eta^4 ...
      + kap*exp(-(x-1.75).^2/eta^2) .* (-2/eta^2);
true4 = kap*exp(-(x-1.75).^2/eta^2) .* (12*(x-1.75)   )/eta^4 ...
      + kap*exp(-(x-1.75).^2/eta^2) .* (-8*(x-1.75).^3)/eta^6;
xx = 7/4; % peak
t1xx  = kap*exp(-(xx-1.75).^2/eta^2);
t2xx  = kap*exp(-(xx-1.75).^2/eta^2) .* (-2*(xx-1.75))/eta^2;
t3xx  = kap*exp(-(xx-1.75).^2/eta^2) .* ( 4*(xx-1.75).^2)/eta^4 ...
      + kap*exp(-(xx-1.75).^2/eta^2) .* (-2/eta^2);
t4xx  = kap*exp(-(xx-1.75).^2/eta^2) .* (12*(xx-1.75)   )/eta^4 ...
      + kap*exp(-(xx-1.75).^2/eta^2) .* (-8*(xx-1.75).^3)/eta^6;
kap2 =  64*kap + 2*t3xx;
kap1 =-208*kap - 7*t3xx;
kap0 =   4*kap - kap2*xx^2 - kap1*xx;
%kap3 = t4xx/6; % should be 0
%kap2 =(t3xx - 6*kap3*xx)/2;
%kap1 = t2xx - 3*kap3*xx^2 - 2*kap2*xx;
%kap0 = t1xx -   kap3*xx^3 -   kap2*xx^2 - kap1*xx;
%%kap0 = -8*kap3 - 4*kap2 - 2*kap1;
%[t1xx t2xx t3xx t4xx], [kap0 kap1 kap2 kap3]
n = length(x);
for i = 1:n
    if (ax <= x(i)) & (x(i) <= 0)
    elseif (0 < x(i)) & (x(i) < 7/4)
	xx = x(i);
	true1(i) = kap*exp(-(xx-1.75).^2/eta^2);
	true2(i) = kap*exp(-(xx-1.75).^2/eta^2) .* (-2*(xx-1.75))/eta^2;
	true3(i) = kap*exp(-(xx-1.75).^2/eta^2) .* ( 4*(xx-1.75).^2)/eta^4 ...
         	 + kap*exp(-(xx-1.75).^2/eta^2) .* (-2/eta^2);
        true4(i) = kap*exp(-(xx-1.75).^2/eta^2) .* (12*(xx-1.75)   )/eta^4 ...
                 + kap*exp(-(xx-1.75).^2/eta^2) .* (-8*(xx-1.75).^3)/eta^6;
    elseif (7/4 <= x(i)) & (x(i) < 2)
        true1(i) = (2-x(i))*(  kap2*x(i)^2 + kap1*x(i) + kap0);
        true2(i) = (2-x(i))*(2*kap2*x(i)   + kap1) ...
                 - (  kap2*x(i)^2 + kap1*x(i) + kap0);
        true3(i) = (2-x(i))*(2*kap2) - 2*(2*kap2*x(i)   + kap1);
        true4(i) = -6*kap2;
        %true1(i) =   kap3*x(i)^3 +   kap2*x(i)^2 + kap1*x(i) + kap0;
        %true2(i) = 3*kap3*x(i)^2 + 2*kap2*x(i)   + kap1;
        %true3(i) = 6*kap3*x(i)   + 2*kap2;
        %true4(i) = 6*kap3;
    elseif (2 <= x(i)) & (x(i) <= bx)
        true1(i) = 0;
        true2(i) = 0;
        true3(i) = 0;
        true4(i) = 0;
    end
end
% incomplete
true5 = zz;
true6 = zz;

case {-103.31}
    Uname ='initial & boundary data only - piece. quadr. with DD at 0, 2 + exp';
true1 = kap*exp(-(x-1.75).^2/eta^2);
true2 = kap*exp(-(x-1.75).^2/eta^2) .* (-2*(x-1.75))/eta^2;
true3 = kap*exp(-(x-1.75).^2/eta^2) .* ( 4*(x-1.75).^2)/eta^4 ...
      + kap*exp(-(x-1.75).^2/eta^2) .* (-2/eta^2);
n = length(x);
for i = 1:n
    if (ax <= x(i)) & (x(i) <= 0)
    elseif (0 < x(i)) & (x(i) < 2)
        true1(i) = true1(i) -(x(i)-1)^2 + 1;
        true2(i) = true2(i) -2*(x(i)-1);
        true3(i) = true3(i) -2;
    elseif (2 <= x(i)) & (x(i) <= bx)
    end
end
true5 = zz;
true6 = zz;

case {-103.32}
    Uname ='initial & boundary data only - mix of exp and cub with DD at 0, 2';
% matching up to including second derivative at 7/4
% discontinuous first derivative at 2, fully conditnuous elsewhere
true1 = zz;
true2 = zz;
true3 = zz;
true4 = zz;
xx = 7/4; % peak
lam = -32; % second derivative at xx
lam2 = lam/2/xx + kap/xx^3;
lam1 =-kap/xx^2 - 2*lam2*xx;
lam0 = kap/xx - xx^2*lam2 - xx*lam1;
n = length(x);
for i = 1:n
    if (ax <= x(i)) & (x(i) <= 0)
    elseif (0 < x(i)) & (x(i) < 7/4)
        xx = x(i);
        true1(i) = xx.*(  lam2*xx^2 + lam1*xx + lam0);
        true2(i) = xx.*(2*lam2*xx   + lam1) + lam2*xx^2 + lam1*xx + lam0;
        true3(i) = xx*2*lam2 + 2*(2*lam2*xx   + lam1);
        true4(i) = 6*lam2;
    elseif (7/4 <= x(i)) & (x(i) < 2)
        true1(i) = (2-x(i))*(  kap2*x(i)^2 + kap1*x(i) + kap0);
        true2(i) = (2-x(i))*(2*kap2*x(i)   + kap1) ...
                 - (  kap2*x(i)^2 + kap1*x(i) + kap0);
        true3(i) = (2-x(i))*(2*kap2) - 2*(2*kap2*x(i)   + kap1);
        true4(i) = -6*kap2;
    elseif (2 <= x(i)) & (x(i) <= bx)
    end
end
% incomplete
true5 = zz;
true6 = zz;

case {-103.4}
    Uname ='initial and boundary data only - piece. quart. with DD at -1, +1';
% dependence on y variable does not make a difference as the BCs are Neumann
true1 = zz;
true2 = zz;
true3 = zz;
n = length(x);
for i = 1:n
    if (ax <= x(i)) & (x(i) <= -1)
    elseif (-1 < x(i)) & (x(i) < 1)
        true1(i) = (-x(i)^4 + 1)*exp(-yy);
        true2(i) = -4*x(i)^3 * exp(-yy);
        true3(i) = -6*x(i)^2 * exp(-yy);
        true4(i) = -12*x(i) * exp(-yy);
        true5(i) = -12 * exp(-yy);
    elseif (1 <= x(i)) & (x(i) <= bx)
    end
end
true4 = zz;
true5 = zz;
true6 = zz;

case {102}
    Uname ='sin(x).*y';
    true1 = sin(x).*y;       true2 = cos(x).*y;       true3=-sin(x).*y;
    true4 = sin(x);          true5 = cos(x);          true6= zz;
    D4X = true1;
case {101}
    Uname ='cos(x)';
    true1 = cos(x);          true2 =-sin(x);          true3 =-cos(x);
    true4 = zz;              true5 = zz;              true6 = zz;
    D4X = true1;
case {100}
    Uname = 'sin(x)';
    true1 = sin(x);          true2 = cos(x);          true3 =-sin(x);
    true4 = zz;              true5 = zz;              true6 = zz;
    D4X = true1;

case {55}
% uxxxx, uyyyy, uxxxxyyyy = 0 on [0, 1] x [0, 1]
Uname = '(x.^6/360-x.^5/120).*(y.^6/360-y.^5/120) -- u4x = u4y = u4xy = 0';
tx1 = x.^6/360-x.^5/120;   tx2 = x.^5/60-x.^4/24;     tx3 = x.^4/12-x.^3/6;
ty1 = y.^6/360-y.^5/120;   ty2 = y.^5/60-y.^4/24;     ty3 = y.^4/12-y.^3/6;
true1 = tx1.*ty1;          true2 = tx2.*ty1;          true3 = tx3.*ty1;
true4 = tx1.*ty2;          true5 = tx2.*ty2;          true6 = tx1.*ty3;
D4X = (x.^2 - x).*ty1;

% uxxxx, uxxxxyyyy = 0 on [0, 1] x [0, 1]
%true1 = tx1*sin(y);       true2 = tx2*sin(y);         true3 = tx3*sin(y);
%true4 = tx1*cos(y);       true5 = tx2*cos(y);         true6 =-tx1*sin(y);
% uyyyy, uxxxxyyyy = 0 on [0, 1] x [0, 1]
%true1 = sin(x)*ty1;       true2 = cos(x)*ty1;         true3 =-sin(x)*ty1;
%true4 = sin(x)*ty2;       true5 = cos(x)*ty2;         true6 = sin(x)*ty3;
% uxxxx, uyy, uxxxxyy, uxxxxyyyy = 0 on [0, 1] x [0, 1]
%ty1 = (y.^4/12-y.^3/6);   ty2 = y.^3/3-y.^2/2;        ty3 = y.^2-y;
%true1 = tx1*ty1;          true2 = tx2*ty1;            true3 = tx3*ty1;
%true4 = tx1*ty2;          true5 = tx2*ty2;            true6 = tx1*ty3;
%
%tx1 = (x.^4/12-x.^3/6);   tx2 = x.^3/3-x.^2/2;        tx3 = x.^2-x;
%ty1 = (y.^4/12-y.^3/6);   ty2 = y.^3/3-y.^2/2;        ty3 = y.^2-y;
%true1 = tx1*ty1;          true2 = tx2*ty1;            true3 = tx3*ty1;
%true4 = tx1*ty2;          true5 = tx2*ty2;            true6 = tx1*ty3;

case {49}
    Uname ='x.^5 .*   y.^5';
    true1 = x.^5 .*   y.^5;  true2 = 5*x.^4 .*   y.^5; true3 =20*x.^3 .*   y.^5;
    true4 = x.^5 .*5.*y.^4;  true5 = 5*x.^4 .*5.*y.^4; true6 =   x.^5 .*20.*y.^3;
case {42}
    Uname ='y.^5';
    true1 = y.^5;            true2 = zz;              true3 = zz;
    true4 = 5*y.^4;          true5 = zz;              true6 = 20*y.^3;
case {39}
    Uname ='x.^4 .*   y.^4';
    true1 = x.^4 .*   y.^4;  true2 = 4*x.^3 .*   y.^4; true3 =12*x.*x .*   y.^4;
    true4 = x.^4 .*4.*y.^3;  true5 = 4*x.^3 .*4.*y.^3; true6 =   x.^4 .*12.*y.*y;
case {32}
    Uname ='y.^4';
    true1 = y.^4;            true2 = zz;              true3 = zz;
    true4 = 4*y.^3;          true5 = zz;              true6 = 12*y.^2;
case {31}
    Uname ='x.^4';
    true1 = x.^4;            true2 = 4*x.^3;          true3 =12*x.*x;
    true4 = zz;              true5 = zz;              true6 = zz;
case {29}
    Uname ='x.^3 .*   y.^3';
    true1 = x.^3 .*   y.^3;  true2 = 3*x.*x .*   y.^3; true3 = 6*x    .*   y.^3;
    true4 = x.^3 .*3.*y.*y;  true5 = 3*x.*x .*3.*y.*y; true6 =   x.^3 .*6.*y;
case {28}
    Uname ='(x.^3 - x.*x)   .*(y.^3 - y.*y)';
    true1 = (x.^3 - x.*x)   .*(y.^3 - y.*y);
    true2 = (3.*x.*x - 2.*x).*(y.^3 - y.*y);
    true3 = (6.*x - 2)      .*(y.^3 - y.*y);
    true4 = (x.^3 - x.*x)   .*(3.*y.*y - 2.*y);
    true5 = (3.*x.*x - 2.*x).*(3.*y.*y - 2.*y);
    true6 = (x.^3 - x.*x)   .*(6.*y - 2);
case {27}
    %ux = 0, uy = 0
    Uname ='(x.^3/3 - x.*x/2).*(y.^3/3 - y.*y/2)';
    true1 = (x.^3/3 - x.*x/2).*(y.^3/3 - y.*y/2);
    true2 = (x.*x - x)       .*(y.^3/3 - y.*y/2);
    true3 = (2*x-1)          .*(y.^3/3 - y.*y/2);
    true4 = (x.^3/3 - x.*x/2).*(y.*y - y);
    true5 = (x.*x - x)       .*(y.*y - y);
    true6 = (x.^3/3 - x.*x/2).*(2*y-1);
case {26}
    Uname ='(x.^3 - x)  .*(y.^3 - y)';
    true1 = (x.^3 - x)  .*(y.^3 - y);
    true2 = (3*x.*x - 1).*(y.^3 - y);
    true3 = 6*x         .*(y.^3 - y);
    true4 = (x.^3 - x)  .*(3.*y.*y - 1);
    true5 = (3*x.*x - 1).*(3*y.*y - 1);
    true6 = (x.^3 - x)  .*6*y;
case {25}
    Uname ='x.*x.*x.*y.*y';
    true1 = x.*x.*x.*y.*y;   true2 = 3*x.*x.*y.*y;    true3 = 6*x.*y.*y;
    true4 = x.*x.*x*2.*y;    true5 = 3*x.*x.*2.*y;    true6 = x.*x.*x*2;
case {22}
    Uname ='y.*y.*y';
    true1 = y.*y.*y.*oo;     true2 = zz;              true3 = zz;
    true4 = 3*y.*y.*oo;      true5 = zz;              true6 = 6*y.*oo;
case {21}
    Uname ='x.*x.*x';
    true1 = x.*x.*x;         true2 = 3*x.*x;          true3 = 6*x;
    true4 = zz;              true5 = zz;              true6 = zz;
case {19}
    Uname ='x.*x.*y.*y';
    true1 = x.*x.*y.*y;      true2 = 2*x.*y.*y;       true3 = 2*y.*y;
    true4 = x.*x.*2.*y;      true5 = 2*x.*2.*y;       true6 = 2*x.*x;
case {18}
    Uname ='(x.*x - x).*(y.*y - y)';
    true1 = (x.*x - x).*(y.*y - y);
    true2 = (2*x - 1) .*(y.*y - y);
    true3 = 2          *(y.*y - y);
    true4 = (x.*x - x).*(2*y - 1);
    true5 = (2*x - 1) .*(2*y - 1);
    true6 = 2          *(x.*x - x);
% u = 0 on x = 0, x = 1, y = 0, y = 1
%true1 = (x*x - x)*(y*y - y); true2 = (2*x - 1)*(y*y - y); true3 = 2*(y*y - y);
%true4 = (x*x - x)*(2*y - 1); true5 = (2*x - 1)*(2*y - 1); true6 = 2*(x*x - x);
% u = 0 on x = 0, x = 1, uy = 0 (on y = 0, y = 1)
%true1 =-(x*x - x)/2; true2 =-(2*x - 1)/2; true3 =-1;
%true4 = 0; true5 = 0; true6 = 0;
% u = 0 on y = 0, y = 1, ux = 0 (on x = 0, x = 1)
%true1 =-(y*y - y)/2; true2 = 0; true3 = 0;
%true4 =-(2*y - 1)/2; true5 = 0; true6 =-1;
case {17}
    Uname ='x.*x + y.*y';
    true1 = x.*x + y.*y;     true2 = 2*x;             true3 = 2*oo;
    true4 = 2.*y.*oo;        true5 = zz;              true6 = 2*oo;
case {16}
    Uname ='x.*y.*y';
    true1 = x.*y.*y;         true2 = y.*y.*oo;        true3 = zz;
    true4 = 2*x.*y;          true5 = 2*y.*oo;         true6 = 2*x;
case {15}
    Uname ='x + y.*y';
    true1 = x + y.*y;        true2 = oo;              true3 = zz;
    true4 = 2.*y.*oo;        true5 = zz;              true6 = 2*oo;
case {14}
    Uname ='y.*y';
    true1 = y.*y.*oo;        true2 = zz;              true3 = zz;
    true4 = 2.*y.*oo;        true5 = zz;              true6 = 2*oo;
case {13}
    Uname ='x.*x.*y';
    true1 = x.*x.*y;         true2 = 2*x.*y;          true3 = 2*oo.*y;
    true4 = x.*x;            true5 = 2*x;             true6 = zz;
case {12}
    Uname ='x.*x + y';
    true1 = x.*x + y;        true2 = 2*x;             true3 = 2*oo;
    true4 = oo;              true5 = zz;              true6 = zz;
case {11}
    Uname ='x.*x';
    true1 = x.*x;            true2 = 2*x;             true3 = 2*oo;
    true4 = zz;              true5 = zz;              true6 = zz;
case {10}
    Uname ='x.*(1-x)';
    true1 = x.*(1-x);        true2 = 1-2*x;           true3 = -2*oo;
    true4 = zz;              true5 = zz;              true6 = zz;
case {9}
    Uname ='x.*y';
    true1 = x.*y;            true2 = y.*oo;           true3 = zz;
    true4 = x;               true5 = oo;              true6 = zz;
case {8}
    Uname ='x.*(1-y)';
    true1 = x.*(1-y);        true2 = 1-y;             true3 = zz;
    true4 =-x;               true5 =-oo;              true6 = zz;
case {5}
    Uname ='x + y';
    true1 = x + y;           true2 = oo;              true3 = zz;
    true4 = oo;              true5 = zz;              true6 = zz;
case {4}
    Uname ='1-y';
    true1 = 1-y;             true2 = zz;              true3 = zz;
    true4 =-oo;              true5 = zz;              true6 = zz;
case {3}
    Uname ='y';
    true1 = y.*oo;           true2 = zz;              true3 = zz;
    true4 = oo;              true5 = zz;              true6 = zz;
case {2}
    Uname ='-x';
    true1 = -x;              true2 =-oo;              true3 = zz;
    true4 =  zz;             true5 = zz;              true6 = zz;
case {1}
    Uname ='x';
    true1 = x;               true2 = oo;              true3 = zz;
    true4 = zz;              true5 = zz;              true6 = zz;
case {0}
    Uname ='1';
    true1 = oo;              true2 = zz;              true3 = zz;
    true4 = zz;              true5 = zz;              true6 = zz;
case {-1}
    Uname ='0';
    true1 = zz;              true2 = zz;              true3 = zz;
    true4 = zz;              true5 = zz;              true6 = zz;
case {-2}
    Uname ='initial and boundary data only: 1 for x <= 0, -1 for x > 0';
    true1 = zz;              true2 = zz;              true3 = zz;
    true4 = zz;              true5 = zz;              true6 = zz;
    true1 = [x <= 0]*2 - 1;
otherwise
    error(['truevd2: no such function ' num2str(Uno)])
end
