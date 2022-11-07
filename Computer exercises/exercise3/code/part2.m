clear all, close all, clc

k = 401;
Cp = 385;
rho = 8963;
a = 1/4;
tpL2 = (a*rho*Cp)/(k)
1/tpL2

n = 10;
xlim = 1;
dx = (xlim-0)/(n+1);
x = dx:dx:xlim;

taulim = 2;
dt = (1/10)*dx^2;
tau = dt:dt:taulim;

A = diag((-2)*ones(1,n), 0) + diag(ones(1,n-1), 1) + diag(ones(1,n-1), -1);
A(n, n+1) = 1;
A(n+1, n) = 2;
A(n+1, n+1) = -2;
A = (a/(dx^2))*A;
A = sparse(A);

b = zeros(n+1,1);

g = zeros(n+1,1);


%% Explicit euler
close all, clc

U = [g];
u = U;
u0 = alph(0);

for t = tau
    b(1) = alph(t)*a/(dx^2);
    u0 = [u0 alph(t)];
    u = u + dt*(A*u + b);
    U = [U u];
end

U = [u0; U];
tau = [0 tau];
x = [0 x];
mesh(U)
title('3D-plot of temperature as function of time and x')
xlabel('time')
ylabel('x')
zlabel('temperature')

dx
dt
dt/(dx^2)

figure
plot(x, U(:,round(1/dt)))
title('Temperature for fixed time \tau = 1')
xlabel('x')
ylabel('temperature')

figure
plot(tau, U(1,:))
hold on
plot(tau, U(end,:))
title('Temperature as function of time at endpoints of rod')
legend('u(0,\tau)', 'u(1,\tau)')
xlabel('time')
ylabel('temperature')


%% ODE23
close all, clc

tspan = [0 2];
y0 = zeros(n+1,1);
bfunc = @(t) [alph(t)*a/(dx^2); zeros(n,1)];
odefunc = @(t,x) A*x+bfunc(t);
opts = odeset('Jacobian',A);

n

tic
[t,y] = ode23(odefunc, tspan, y0);
toc
Timesteps = length(t)
plot(t,y)

tic
[ts,ys] = ode23s(odefunc, tspan, y0);
toc
Timesteps = length(ts)
plot(t,y)


tic
[ts,ys] = ode23s(odefunc, tspan, y0, opts);
toc
Timesteps = length(ts)
plot(t,y)








