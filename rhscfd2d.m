function [ rhs, coefs ] = rhscfd2d( gridx, gridy )
%RHSCFD return the right hand side vector of solvinh PDE
%   the RHS vector is obtained from the combination of RHS function g(x,y)
%   and the boundary conditions
%   Dirichlet conditions, boundary conditions stored in bc2.m
%   boundaries are {x = ax, x = bx, y = ay, y = by}

n = length(gridx) - 1; % number of grid on x
m = length(gridy) - 1; % number of grid on y
numeq = (n - 1)*(m - 1);	% assume Dirichlet conditions
hx = gridx(2:end) - gridx(1:end-1);
hy = gridy(2:end) - gridy(1:end-1);

%boundary values from Dirichlet BSc
for i = 1:(n-1)
    bv1(i) = bc2(gridx(end), gridy(i+1),1);
    bv3(i) = bc2(gridx(1), gridy(i+1),3);
end
for i = 1:(m-1)
    bv2(i) = bc2(gridx(i+1), gridy(1),2);
    bv4(i) = bc2(gridx(i+1), gridy(end),4);
end
bv5 = bc2(gridx(end), gridy(1), 5); bv6 = bc2(gridx(1), gridy(1), 6);
bv7 = bc2(gridx(1), gridy(end), 7); bv8 = bc2(gridx(end), gridy(end), 8);

%RHS , COEFs and Modification of RHS from BCs
counter = 1;
for i = 1:(n-1)
    for j = 1:(m-1)
        px = gridx(i+1); py = gridy(j+1);
        [rhs(counter, 1), coefs(counter, 1), coefs(counter, 2), ...
            coefs(counter, 3), coefs(counter, 4), coefs(counter, 5),...
            coefs(counter, 6)] = PDEcoefs(px, py);

        %modification of RHS from the boundary condition
        
        %u_xx and u_x
        if i == 1
            rhs(counter) = rhs(counter)-bv3(j)*...
                (2*coefs(j,3) - coefs(j,2)*hx(2))/(hx(1)*(hx(1)+hx(2)));
        elseif i == n-1
            rhs(counter) = rhs(counter)-bv1(j)*...
                (2*coefs(counter,3) + coefs(counter,2)*hx(end-1))/(hx(end)*(hx(end)+hx(end-1)));
        end
        
        %u_yy and u_y
        if j == 1
            rhs(counter) = rhs(counter)-bv2(i)*...
                (2*coefs(counter,6) - coefs(counter,4)*hy(2))/(hy(1)*(hy(1)+hy(2)));
        elseif j == m-1
            rhs(counter) = rhs(counter)-bv4(i)*...
                (2*coefs(counter,6) + coefs(counter,4)*hy(end-1))/(hy(end)*(hy(end)+hy(end-1)));
        end
        
        %u_xy
        if i == 1 && j == 1
            rhs(counter) = rhs(counter) - coefs(counter, 5)*...
                (bv6-bv3(2)-bv2(2))/(hx(1)+hx(2))/(hy(1)+hy(2));
        elseif i == 1 && j == (m-1)
            rhs(counter) = rhs(counter) - coefs(counter, 5)*...
                (bv3(end-1)+bv4(2)-bv7)/(hx(1)+hx(2))/(hy(end)+hy(end-1));
        elseif i == (n-1) && j == 1
            rhs(counter) = rhs(counter) - coefs(counter, 5)*...
                (bv2(end-1)+bv1(2)-bv5)/(hy(1)+hy(2))/(hx(end)+hx(end-1));
        elseif i == (n-1) && j == (m-1)
            rhs(counter) = rhs(counter) - coefs(counter, 5)*...
                (bv8-bv4(end-1)-bv1(end-1))/(hx(end)+hx(end-1))/(hy(end)+hy(end-1));
        elseif i == 1 && j ~= 1 && j ~= (m-1)
            rhs(counter) = rhs(counter) - coefs(counter, 5)*...
                (bv3(j-1)-bv3(j+1))/(hx(1)+hx(2))/(hy(j)+hy(j+1));
        elseif i == (n-1) && j ~= 1 && j ~= (m-1)
            rhs(counter) = rhs(counter) - coefs(counter, 5)*...
                (-bv1(j-1)+bv1(j+1))/(hx(end)+hx(end-1))/(hy(j)+hy(j+1));
        elseif j == 1 && i ~= 1 && i ~= (n-1)
            rhs(counter) = rhs(counter) - coefs(counter, 5)*...
                (bv2(i-1)-bv2(i+1))/(hx(i)+hx(i+1))/(hy(1)+hy(2));
        elseif j == (m-1) && i ~= 1 && i ~= (n-1)
            rhs(counter) = rhs(counter) - coefs(counter, 5)*...
                (-bv4(i-1)+bv4(i+1))/(hx(i)+hx(i+1))/(hy(end-1)+hy(end));
        end
       
        counter = counter + 1;
    end
end
end
