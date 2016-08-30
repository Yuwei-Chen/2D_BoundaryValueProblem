%Test function setup
format compact
global Uno Uname;
Uno = 116; Uname = '';

%Initial setup
ax = 0; bx = 1; % domain on x
ay = 0; by = 1; % domain on y
ntimes = 4;	% number of grid sizes to run
errg = zeros(1, ntimes);

%solving PDEs for different size of grid for convergence check
for nn = 1:ntimes
    n = 2^(nn+2); ngrid = n+1; nint(nn) = n;
    neq = (n-1)^2; numeq = neq;	% (n-1)*(m-1) for D
    
    %set up of non-uniform grid pn x and y
    hx_u = (bx-ax)/n;
    gridx_u = ax + hx_u*[0:n];
    gridx = gridx_u.^2;
    hy_u = (by-ay)/n;
    gridy_u = ax + hx_u*[0:n];
    gridy = gridy_u.^2;

    %solving PDE
    [rhs, coefs] = rhscfd2d(gridx, gridy);
    [ A ] = cfdmat2d(gridx, gridy, coefs);
    uvct = A\rhs;
    
    %Inf-norm abs error (approximation and true value)
    errg = errorfd(ngrid, gridx, gridy, n, uvct, nn, errg);
end

%OUTPUT
[udummy] = truevd2(ax,ay);
disp(['U = ' Uname ' = {' num2str(Uno) '}'])
disp(['domain [' num2str(ax)  ', ' num2str(bx) '] X [' num2str(ax)  ', ' num2str(bx) ']']);

nint
format short e
disp('error on grid points')
errg
format short
disp('order of convergence')
errg(:, :) = max(errg(:, :),  0.222044604925e-15);
LogNintRatio = log(nint(1, 2:ntimes)./nint(1, 1:ntimes-1));
LogNintRatioMat = repmat(LogNintRatio, size(errg, 1), 1);
if ntimes > 1
    convg = log(errg(:, 1:ntimes-1)./errg(:, 2:ntimes))./LogNintRatioMat
end