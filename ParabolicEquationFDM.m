clc 
clear

L = 100;
T0 = 100;
a = 1;


% Euler Method
tmax = 1000;
dx = 0.5;
dt = 0.1;
r = dt/(a*dx^2);

x = 0:dx:L;
t = 0:dt:tmax;

u_euler = zeros(length(x), length(t));
u_euler(:, 1) = T0*(sin(pi*x/L).^2);


for j = 2:length(t)
    for i = 2:length(x)-1
        u_euler(i, j) = r*u_euler(i+1, j-1) + (1-2*r)*u_euler(i, j-1) + r*u_euler(i-1, j-1);
    end
end

X = linspace(0, L, length(x));
T = linspace(0, tmax, length(t));
[X,T] = meshgrid(X, T);

figure(1);
contourf(X, T, u_euler');
colorbar;
title("Parabolic Equation Solution with Euler Method");


% Crank-Nicolson
tmax = 1000;
dx = 0.4;
dt = 0.1;
r = dt/(a*dx^2);
lambda = 0.4;

x = 0:dx:L;
t = 0:dt:tmax;
M = length(x);
N = length(t);

u_CN = zeros(N, M);
u_CN(1, :) = T0*(sin(pi*x/L).^2);

A = (1+2*(1-lambda)*r)*diag(ones(M-1,1))-(1-lambda)*r*diag(ones(M-2,1),1)...
    -(1-lambda)*r*diag(ones(M-2,1),-1);
B = zeros(M-1, 1);

for n = 2 : N 
    for m = 2 : M-1 
         B(m-1) = lambda*r*u_CN(n-1,m+1)+(1-2*lambda*r)*u_CN(n-1,m)...
             +lambda*r*u_CN(n-1,m-1);
    end

    sol = A\B; 
    u_CN(n,2:M) = sol'; 
end

X = linspace(0, L, length(x));
T = linspace(0, tmax, length(t));
[X,T] = meshgrid(X, T);

figure(2);
contourf(X, T, u_CN);
colorbar;
title("Parabolic Equation Solution with Crank-Nicolson Method")