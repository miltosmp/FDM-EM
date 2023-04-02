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
    for i = 1:length(x)
        if (i == 1) || (i == length(x))
            u_euler(i, j) = 0;
        else
            u_euler(i, j) = r*u_euler(i+1, j-1) + (1-2*r)*u_euler(i, j-1) + r*u_euler(i-1, j-1);
        end
    end
end

X = linspace(0, L, length(x));
T = linspace(0, tmax, length(t));
[X,T] = meshgrid(X, T);

figure(1);
contourf(X, T, u_euler');
colorbar;


