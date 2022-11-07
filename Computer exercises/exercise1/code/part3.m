clear all, close all, clc

r1 = 0.04;
r2 = 1e4;
r3 = 3e7;

dx = @(x1, x2, x3) [((-r1)*x1) + (r2*x2*x3); (r1*x1) - (r2*x2*x3) - (r3*(x2^2)); (r3*(x2^2))];

%% RK3


N = [125 250 500 1000 2000];
T = 1;
c = 1;

for N = N
    x = [1; 0; 0];
    h = T/N;
    t = 0:h:T;


    x = rk3p3(x, t, h, dx);

    figure
    loglog(t, x(1,:))
    hold on
    loglog(t, x(2,:))
    hold on
    loglog(t, x(3,:))
    legend('x1','x2','x3')
    title("N = " + N)
    xlabel('log(t)')
    ylabel('log(x)')
    c = c+1;
end



%% ODE23
clc

T = 1;
tspan = [0 T];
x = [1; 0; 0];

rel = [1e-3 1e-4 1e-5 1e-11];

c = 1;
for rel = rel
    
    abs = rel/1000;
    options = odeset('RelTol', rel, 'AbsTol', abs);

    [t,xode] = ode23('dxODE', tspan, x, options);
    

    
%     subplot(3,2,c);
%     loglog(t, xode(:,1))
%     hold on
%     loglog(t, xode(:,2))
%     hold on
%     loglog(t, xode(:,3))
%     legend('x1','x2','x3')
%     title(rel)
%     c = c+1;

    if c == 3 | c == 4
        hdiff = diff(t);
        plot(t(1:end-1), hdiff)
        hold on
        length(t)
    end
    
    c = c+1;
end
legend('rel=1e-5', 'rel=1e-11')
title('Time difference (h) as function of time (t), ODE23')
xlabel('t')
ylabel('h')



%% ODES23


T = 1000;
tspan = [0 T];
x = [1; 0; 0];

rel = [1e-3 1e-4 1e-5 1e-11];

c = 1;
for rel = rel
    
    abs = rel/1000;
    options = odeset('RelTol', rel, 'AbsTol', abs);

    [t,xode] = ode23s('dxODE', tspan, x, options);
    

    
%     subplot(3,2,c);
%     loglog(t, xode(:,1))
%     hold on
%     loglog(t, xode(:,2))
%     hold on
%     loglog(t, xode(:,3))
%     legend('x1','x2','x3')
%     title(rel)
%     c = c+1;

    if c == 3 | c == 4
        hdiff = diff(t);
        plot(t(1:end-1), hdiff)
        hold on
        length(t)
    end
    
    c = c+1;
end
legend('rel=1e-5', 'rel=1e-11', 'Location', 'northwest')
title('Time difference (h) as function of time (t), ODES23')
xlabel('t')
ylabel('h')




%% Functions


function x = rk3p3(x,t,h,dx)
    for n = 1:(length(t)-1)

        k1 = dx(x(1,n), x(2,n), x(3,n));
        k2 = dx(x(1,n) + h*k1(1), x(2,n) + h*k1(2), x(3,n) + h*k1(3));
        k3 = dx(x(1,n) +((h*k1(1))/4) + ((h*k2(1))/4), x(2,n) +((h*k1(2))/4) + ((h*k2(2))/4), x(3,n) +((h*k1(3))/4) + ((h*k2(3))/4));

        x(:,n+1) = x(:,n) + (h/6)*(k1 + k2 + 4*k3);

    end
end

