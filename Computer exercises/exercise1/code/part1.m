clear all, close all, clc

alpha = 0.1;
a = 1/2*[1,1,sqrt(2)]';

h = 2.1;
T = 100;

t = 0:h:T;
m = [1,0,0]';

dmdt = @(m)(cross(a, m) + alpha .* (cross(a, cross(a, m))));


%% 2D-plot

m = rk3(m,t,h,dmdt);
plot(t, m)
title('Components of m as a function of t (h = 2.1, T = 40)')
xlabel('time')
ylabel('value of components of m')

%% 3D-plot

absm = vecnorm(m);
absm = absm.^(-1);
m3d = m*diag(absm);
m3d(isnan(m3d)) = 0;

plot3([m3d(1,:),0], [m3d(2,:),0], [m3d(3,:),0])
hold on
plot3([0,a(1)],[0,a(2)],[0,a(3)])

axis equal
legend('m/|m|','a')
title('Trajectory of m/|m|')


%% Order of accuracy

diff = [];
N = [20,40,80,160,320];

for n = N
    m1 = [1,0,0]';
    h1 = T/n;
    t1 = 0:h1:T;
    m1 = rk3(m1,t1,h1,dmdt);
    m1 = m1(:,end);
    
    m2 = [1,0,0]';
    h2 = T/(2*n);
    t2 = 0:h2:T;
    m2 = rk3(m2,t2,h2,dmdt);
    m2 = m2(:,end);
    
    diff = [diff, vecnorm(m1-m2)];
end

hval = (zeros(1,length(N)) + T)./N;


loglog(hval, diff)
hold on
loglog(hval, hval.^3)
hold on
loglog(hval, diff, 'o')

legend('|mN - m2N|', 'h^3')
title('Order of accuracy in loglog')
xlabel('stepsize')
ylabel('value')



%% Step size stability limit
clc

A = (1/2)*[(-3/20), ((1/20) - sqrt(2)), ((sqrt(2)/20) + 1);
    ((1/20) + sqrt(2)), (-3/20), ((sqrt(2)/20) - 1);
    ((sqrt(2)/20) -1), ((sqrt(2)/20) + 1), (-1/10)];

e = eig(A);
e = e(1);
hstab = 0.05;
z = e*hstab;
check = norm(1 + z + (z^2/2) + (z^3/6));

while check < 1
    hstab = hstab + 0.0001;
    z = e*hstab;
    check = norm(1 + z + (z^2/2) + (z^3/6));
end

hstab



%% Functions


function m = rk3(m,t,h,dmdt)
    for n = 1:(length(t)-1)

        k1 = dmdt(m(:,n));
        k2 = dmdt(m(:,n) + h*k1);
        k3 = dmdt(m(:,n) +((h*k1)/4) + ((h*k2)/4));

        m(:,n+1) = m(:,n) + (h/6)*(k1 + k2 + 4*k3);
        
    end
end

