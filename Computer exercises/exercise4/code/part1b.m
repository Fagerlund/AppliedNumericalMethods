clear all, close all, clc


tau = 2.2;
D = 4.5;
a = 1;
tend = 6;

Nx = 100;
dx = D/(Nx-1);
% dt <= dx
dt = 0.7*dx;
Nt = tend/dt;
lambda = dt/dx;

% x,t - interval
xvec = 0:dx:D;
tvec = 0:dt:tend;

% Boundary conditions
gsin = @(t) sin((2*pi*t)/tau);
gsq = @(t) square((2*pi*t)/tau);

%% b) Sin wave
close all, clc


for lambda = 0.01:0.1:1
    dt = lambda*dx;
    Nt = tend/dt;
    lambda = dt/dx;
    
    % Initial condition: u(x,0) = 0
    % Upwind: uu; Lax-Friedrichs: ulf; Lax-Wendroff: ulw
    uu = zeros(Nx,1);
    ulf = uu;
    ulw = uu;



    for n = 1:Nt
        uu(1,n+1) = gsin(n*dt);
        ulf(1,n+1) = gsin(n*dt);
        ulw(1,n+1) = gsin(n*dt);

        for k = 2:(Nx-1)
            % a > 0, reason for Upwind
            uu(k,n+1) = uu(k,n) - a*lambda*(uu(k,n)-uu(k-1,n));
            ulf(k,n+1) = (1/2)*(ulf(k+1,n)+ulf(k-1,n))-a*(lambda/2)*(ulf(k+1,n)-ulf(k-1,n));
            ulw(k,n+1) = ulw(k,n) - a*(lambda/2)*(ulw(k+1,n)-ulw(k-1,n)) + (a^2)*((lambda^2)/2)*(ulw(k+1,n)-2*ulw(k,n)+ulw(k-1,n));
        end

        uu(k+1,n+1) = 2*uu(k,n+1)-uu(k-1,n+1);
        ulf(k+1,n+1) = 2*ulf(k,n+1)-ulf(k-1,n+1);
        ulw(k+1,n+1) = 2*ulw(k,n+1)-ulw(k-1,n+1);
    end

    figure(1)
    plot(xvec,uu(:,end))
    hold on
    title('Upwind')
    
    figure(2)
    plot(xvec,ulf(:,end))
    hold on
    title('Lax-Friedrich')
    
    figure(3)
    plot(xvec,ulw(:,end))
    hold on
    title('Lax-Wendroff')
    
end


%% b) Square wave
close all, clc


for dt = 0.2:0.2:1
    dt = dt*dx;
    Nt = tend/dt;
    lambda = dt/dx
    
    % Initial condition: u(x,0) = 0
    % Upwind: uu; Lax-Friedrichs: ulf; Lax-Wendroff: ulw
    uu = zeros(Nx,1);
    ulf = uu;
    ulw = uu;

    for n = 1:Nt
        uu(1,n+1) = gsq(n*dt);
        ulf(1,n+1) = gsq(n*dt);
        ulw(1,n+1) = gsq(n*dt);

        for k = 2:(Nx-1)
            % a > 0, reason for Upwind
            uu(k,n+1) = uu(k,n) - a*lambda*(uu(k,n)-uu(k-1,n));
            ulf(k,n+1) = (1/2)*(ulf(k+1,n)+ulf(k-1,n))-a*(lambda/2)*(ulf(k+1,n)-ulf(k-1,n));
            ulw(k,n+1) = ulw(k,n) - a*(lambda/2)*(ulw(k+1,n)-ulw(k-1,n)) + (a^2)*((lambda^2)/2)*(ulw(k+1,n)-2*ulw(k,n)+ulw(k-1,n));
        end

        uu(k+1,n+1) = 2*uu(k,n+1)-uu(k-1,n+1);
        ulf(k+1,n+1) = 2*ulf(k,n+1)-ulf(k-1,n+1);
        ulw(k+1,n+1) = 2*ulw(k,n+1)-ulw(k-1,n+1);
    end


    figure(1)
    plot(xvec,uu(:,end))
    hold on
    title('Upwind')
    
    figure(2)
    plot(xvec,ulf(:,end))
    hold on
    title('Lax-Friedrich')
    
    figure(3)
    plot(xvec,ulw(:,end))
    hold on
    title('Lax-Wendroff')
end


