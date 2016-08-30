% function errg = errorfd(ngrid, gridx, n, uvct, nn, errg, truef);
function errg = errorfd(ngrid, gridx, gridy, n, uvct, nn, errg, truef)

if (nargin < 8); truef = 'truevd2'; end;
n = length(gridx);
m = length(gridy);
tr=[];
for i = 2:(n-1)
    for j = 2:(m-1)
        [t1] = feval(truef, gridx(i), gridy(j));
        tr = [tr, t1];
    end
end
errg(1, nn) = max(abs(tr-uvct'));