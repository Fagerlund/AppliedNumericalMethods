clear all, close all, clc

T0l = 40;
T5l = 400;

Lx = 5;
Ly = 2;
h = 0.1;
N = (Lx/h)-1;
M = (Ly/h)+1;


%% A - matrix

A = diag([4+zeros(1,N*M)], 0);

A = A + diag(-1+zeros(1,N*M-1), 1);
A = A + diag(-1+zeros(1,N*M-1), -1);

A = A + diag(-1+zeros(1,N*(M-1)), N);
A = A + diag(-1+zeros(1,N*(M-1)), -N);

for i = 1:M-1
    A(i*N,i*N+1) = 0;
    A(i*N+1,i*N) = 0;
end

for n = 1:N
    A(n, N+n) = -2;
    A(N*(M-1)+n, N*(M-2)+n) = -2;
end

A = (1/h^2)*A;


%% gx - matrix

bcx = zeros(N,1);
bcx(1) = T0l;
bcx(end) = T5l;
gx = bcx;
for i = 1:M-1
    gx = [gx; bcx];
end

gx = (1/h^2)*gx;



%% f - matrix

f = 100 + zeros(M*N, 1);
f = f + gx;


%% solve T
clc

T = A\f;
T = reshape(T, N, M);
T = [40 + zeros(1,size(T,2)); T; 400 + zeros(1,size(T,2))];

mesh(T)
title('Mesh-plot of temperature as function of x and y (h=0.1)')
xlabel('y')
ylabel('x')
zlabel('T - temp')


T((3/h)+1, (1/h)+1)

% (3,1) = 556


%% c)
clc, close all


ff = @(x,y) 6000*exp(-5*((x-1)^2) - 10*((y-1.5)^2));

fc = zeros(M*N, 1);
for m = 0:M-1
    for n = 1:N
        fc(N*m + n) = ff(n*h,m*h);
    end
end

fc = fc + gx;


T = A\fc;
T = reshape(T, N, M);
T = [40 + zeros(1,size(T,2)); T; 400 + zeros(1,size(T,2))];
T((3/h)+1, (1/h)+1)

% (x,y) = (3,1)
% h = 0.1: 781.7898
% h = 0.05: 782.2751
% h = 0.025: 782.3972

% såg att det går mot 782.42 med mmindre steglängd


T = T';

mesh(T)
set(gca,'YDir','normal')
title('Mesh-plot of temperature as function of x and y (h=0.05)')
xlabel('x')
ylabel('y')
zlabel('T - temp')

figure
imagesc(T)
colorbar
set(gca,'YDir','normal')
title('Imagesc-plot of temperature as function of x and y (h=0.05)')
xlabel('x')
ylabel('y')
zlabel('T - temp')



figure
contour(T)
set(gca,'YDir','normal')
title('Contour-plot of temperature as function of x and y (h=0.05)')
xlabel('x')
ylabel('y')
zlabel('T - temp')





