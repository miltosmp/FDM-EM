clc
clear

a = 2;
b = 1;

M = 100;
N = 50;
dx = a/M;
dy = b/N;

x = linspace(0, a, M);
y = linspace(0, b, N);
[X, Y] = meshgrid(x, y);

n_un = (M-1)*(N-1);

L = ones(n_un-1,1);   
for k = 1 : n_un-1
    if mod(k,M-1) == 0 
        L(k) = 0; 
    end   
end

A = (2/dx^2 + 2/dy^2)*diag(ones(n_un, 1), 0) - (1/dx^2)*diag(L, 1)...
    - (1/dx^2)*diag(L, -1) - (1/dy^2)*diag(ones((M-1)*(N-2), 1), M-1) - ...
    (1/dy^2)*diag(ones((M-1)*(N-2), 1), -M+1);

A = sparse(A);
[V, D] = eigs(A, 10, 'sm');

x = 0:dx:a;
y = 0:dy:b;

for k = 1:10
    U = zeros(N+1, M+1);
    l = 1;
    for j = 1:N-1
        for i = 1:M-1
            U(j+1, i+1) = V(l, k);
            l = l + 1;
        end
    end
    figure(k);
    contourf(x, y, U);
end
