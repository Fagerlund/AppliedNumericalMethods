clear all, close all, clc

u0 = 40;
u5 = 400;

Lx = 5;
Ly = 2;
h = 0.1;

N = (Lx/h)-1;
M = (Ly/h)+1;

xvec = h:h:Lx-h;
yvec = h:h:Ly+h;

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

A = -(1/h^2)*A;
A = sparse(A);


%% gx - matrix and f - matrix

bcx = zeros(N,1);
bcx(1) = u0;
bcx(end) = u5;
gx = bcx;
for i = 1:M-1
    gx = [gx; bcx];
end

gx = (1/h^2)*gx;


ff = @(x,y) 6000*exp(-5*((x-1)^2) - 10*((y-1.5)^2));

fc = [];
for y = yvec
    for x = xvec
        fc = [fc; ff(x,y)];
    end
end

b = fc + gx;



%% Solve u (Crank-Nicolson)

close all, clc

I = eye(N*M);
I = sparse(I);
taulim = 10;
C = 1;
dt = C*(h^2);
tau = dt:dt:taulim;

%u0 (t=0)
u0func = @(x,y) 40+72*x;
u = [];
for y = yvec
    for x = xvec
        u = [u; u0func(x,y)];
    end
end

% Crank-Nicolson
U = u;
for t = tau
    u = (sparse((I-(1/2)*dt.*A)))\( sparse(((I+((1/2)*dt.*A))*u) + dt.*b) );
    
    U = [U u];
end



Ut = U(:,1);
Ut = reshape(Ut, N, M);
Ut = [u0 + zeros(1,size(Ut,2)); Ut; u5 + zeros(1,size(Ut,2))];
mesh(Ut)
title('t = 0')
xlabel('y')
ylabel('x')
zlabel('temp')

figure
Ut = U(:,round(0.25/dt)+1);
Ut = reshape(Ut, N, M);
Ut = [u0 + zeros(1,size(Ut,2)); Ut; u5 + zeros(1,size(Ut,2))];
mesh(Ut)
title('t = 0.25')
xlabel('y')
ylabel('x')
zlabel('temp')

figure
Ut = U(:,round(1/dt)+1);
Ut = reshape(Ut, N, M);
Ut = [u0 + zeros(1,size(Ut,2)); Ut; u5 + zeros(1,size(Ut,2))];
mesh(Ut)
title('t = 1')
xlabel('y')
ylabel('x')
zlabel('temp')

figure
Ut = U(:,round(10/dt)+1);
Ut = reshape(Ut, N, M);
Ut = [u0 + zeros(1,size(Ut,2)); Ut; u5 + zeros(1,size(Ut,2))];
mesh(Ut)
title('t = 10')
xlabel('y')
ylabel('x')
zlabel('temp')

figure
Ut = U(N/h+3/h,:);
tau = [0 tau];
plot(tau, Ut)
hold on
plot(tau, 782.4*ones(length(tau),1))
legend('(3,1) over time', '2c', 'Location', 'Best') 
title('(x,y) = (3,1)')
xlabel('time')
ylabel('temp')






