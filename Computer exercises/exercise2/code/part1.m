

%% a)
clear all, close all, clc

v = 10;
% h = 0.05;
hvec = [0.05 0.025 0.0125 0.00625];
% hvec = [0.003125 0.003125/2 0.003125/4 0.003125/8];


diff = [];

for h = hvec
    [z, T] = solve(h,v);
    [z2, T2] = solve(h/2,v);
    
    diff = [diff, abs(T(end)-T2(end))];
    
    plot(z, T)
    hold on
end

xlabel('z - length')
ylabel('T - temp')
legend('h=0.05', 'h=0.025', 'h=0.0125', 'h=0.00625')
title('T(z) through the pipe with various h')


figure
loglog(hvec, diff)
hold on
loglog(hvec, hvec.^2)
legend('diff', 'h^2')
title('Convergence of solution')
xlabel('h')



%% b)
clear all, close all, clc

v = [1 10 30 100];
% h = 0.00625;
h = 0.03;

for v = v
    
    [z, T] = solve(h,v);
    plot(z, T)
    hold on
end

xlabel('z - length')
ylabel('T - temp')
legend('v=1', 'v=10', 'v=30', 'v=100')
title('T(z) through the pipe with various v')



%% function

function [z, T] = solve(h,v)
    L = 1;
    a = 0.2;
    b = 0.3;
    Q0 = 4000;
    alpha0 = 100;
    Tout = 20;
    T0 = 50;


    alpha = @(v) sqrt(((v^2)/4) + (alpha0^2)) - (v/2);
    
    
    z = h:h:L;


    aj = -1/(h^2) - (v/(2*h));
    bj = 2/(h^2);
    cj = -1/(h^2) + (v/(2*h));


    A = diag([bj+zeros(1,length(z)-1) 0], 0);
    A = A + diag(cj+zeros(1,length(z)-1), 1);
    A = A + diag([aj+zeros(1,length(z)-2) 0], -1);
    A(length(z), length(z)-1) = aj+cj;
    A(length(z), length(z)) = bj - cj*(2*h*alpha(v));            


    q = Q(z(1)) - aj*T0;


    for i = z(2:end-1)
        q = [q; Q(i)];
    end

    q = [q; Q(z(end)) - cj*2*h*alpha(v)*Tout];


    T = A\q;
    T = [T0; T];
    z = [0 z]';

    
end




