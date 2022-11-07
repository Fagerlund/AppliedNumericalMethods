clear all, close all, clc

my = 1/82.45;
r0 = [-my, 0]';
r1 = [1-my, 0]';

N = [0, 1; -1 0];

h = 0.001;
t = 0:h:7;

uvec = @(u1,u2) [u2, (-(1-my)*((u1-r0)./(vecnorm(u1-r0).^3))) - (my*((u1-r1)./(vecnorm(u1-r1).^3))) + (2*N*u2) + u1];


TOL = 0.25;



%% Plots
clc

u1 = [1.2, 0]';
u2 = [0, -1]';


for n = 1:3
    [u1, u2] = RK3(u1, u2, h, n);
end


for n = 4:length(t)-1
    [u1, u2] = adam(u1, u2, h, n);   
end


plot3(u1(1,:), u1(2,:), t)
hold on
plot3(u2(1,:), u2(2,:), t)
hold on
plot3(zeros(1,length(t))+r0(1), zeros(1,length(t))+r0(2), t)
hold on
plot3(zeros(1,length(t))+r1(1), zeros(1,length(t))+r1(2), t)

axis equal
legend('r(t)', 'dr/dt', 'earth', 'moon')
title('3D-plot of position cordinates and velocities over time')
xlabel('x-value')
ylabel('y-value')
zlabel('time')


figure

plot(u1(1,:), u1(2,:))
hold on
plot(u2(1,:), u2(2,:))
hold on
scatter(r0(1), r0(2))
hold on
scatter(r1(1), r1(2))

axis equal
legend('r(t)', 'dr/dt', 'earth', 'moon')
title('2D-plot of position cordinates and velocities (t ? [0, 7])')
xlabel('x-value')
ylabel('y-value')



%% Euler
clc

u1 = [1.2, 0]';
u2 = [0, -1]';

u12 = [1.2, 0]';
u22 = [0, -1]';

boll = true;
n = 1;

while vecnorm(u1(:,end) - u12(:,end)) <= TOL | boll
    boll = false;
    
    [u1, u2] = eul(u1, u2, h, n);
    
    [u12, u22] = eul(u12, u22, h/2, 2*n-1);
    [u12, u22] = eul(u12, u22, h/2, 2*n);
    
    n = n+1;
end

TaccEUL = n * h


% comet(u1(1,:), u1(2,:))
plot(u1(1,:), u1(2,:))
hold on
scatter(r0(1), r0(2))
hold on
scatter(r1(1), r1(2))

legend('r(t)', 'earth', 'moon', 'Location', 'northwest')
title('2D-plot of position cordinates, Euler (t ? [0, Tacc])')
xlabel('x-value')
ylabel('y-value')



%% RK3
clc

u1 = [1.2, 0]';
u2 = [0, -1]';

u12 = [1.2, 0]';
u22 = [0, -1]';

boll = true;
n = 1;

while vecnorm(u1(:,end) - u12(:,end)) <= TOL | boll
    boll = false;

    [u1, u2] = RK3(u1, u2, h, n);
    
    [u12, u22] = RK3(u12, u22, h/2, 2*n-1);
    [u12, u22] = RK3(u12, u22, h/2, 2*n);
    
    
    n = n+1;
end

TaccRK3 = n * h


% comet(u1(1,:), u1(2,:))
plot(u1(1,:), u1(2,:))
hold on
scatter(r0(1), r0(2))
hold on
scatter(r1(1), r1(2))

legend('r(t)', 'earth', 'moon')
title('2D-plot of position cordinates, RK3 (t ? [0, Tacc])')
xlabel('x-value')
ylabel('y-value')



%% Adams
clc

u1 = [1.2, 0]';
u2 = [0, -1]';

u12 = [1.2, 0]';
u22 = [0, -1]';


for n = 1:3
    [u1, u2] = RK3(u1, u2, h, n);
    
    [u12, u22] = RK3(u12, u22, h/2, 2*n-1);
    [u12, u22] = RK3(u12, u22, h/2, 2*n);
end


boll = true;
n = 4;

while vecnorm(u1(:,end) - u12(:,end)) <= TOL | boll
    boll = false;

    [u1, u2] = adam(u1, u2, h, n);
    
    [u12, u22] = adam(u12, u22, h/2, 2*n-1);
    [u12, u22] = adam(u12, u22, h/2, 2*n);
    
    
    n = n+1;
end


% Tacc for adam was like we expected higher than RK and below 20 like it
% said in the description
TaccADAM = n * h


% comet(u1(1,:), u1(2,:))
plot(u1(1,:), u1(2,:))
hold on
scatter(r0(1), r0(2))
hold on
scatter(r1(1), r1(2))

legend('r(t)', 'earth', 'moon')
title('2D-plot of position cordinates, Adams-Bashford (t ? [0, Tacc])')
xlabel('x-value')
ylabel('y-value')



%% ODE23
clc

% uvecODE = @(u1, u2, u3, u4) [u3; u4; (-(1-my)*(([u1; u2]-r0)./(vecnorm([u1; u2]-r0).^3))) - (my*(([u1; u2]-r1)./(vecnorm([u1; u2]-r1).^3))) + (2*N*[u3; u4]) + [u1; u2]];

u1 = [1.2, 0]';
u2 = [0, -1]';

tspan = [0 TaccADAM];
init = [u1(1); u1(2); u2(1); u2(2)];
options = odeset('RelTol', 1e-4);

[t, A] = ode23('uvecODE', tspan, init, options);
    
hdiff = diff(t);

hmax = max(hdiff)
hmin = min(hdiff)
tend = t(end)
tnumb = length(t);
plot(t(1:end-1), hdiff)
title('Time difference (h) as function of time (t)')
xlabel('t')
ylabel('h')

% comet(A(:,1), A(:,2))



