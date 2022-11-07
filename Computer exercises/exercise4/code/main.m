
%% Part 1

clear all, close all, clc


tau = 2.2;
D = 4.5;
a = 1;
tend = 6;

Nx = 100;
dx = D/(Nx-1);
% dt <= dx
dt = 0.4*dx;
Nt = tend/dt;
lambda = dt/dx;

% x,t - interval
xvec = 0:dx:D;
tvec = 0:dt:tend;

% Boundary conditions
gsin = @(t) sin((2*pi*t)/tau);
gsq = @(t) square((2*pi*t)/tau);

%% 1.a) Sin wave
close all, clc

% Initial condition: u(x,0) = 0
% Upwind: uu; Lax-Friedrichs: ulf; Lax-Wendroff: ulw
uu = zeros(Nx,1);
ulf = uu;
ulw = uu;
uuex = uu;


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

dtex = dx;
Ntex = tend/dtex;
lambdaex = dtex/dx;

% exact solution (lambda = 1)
for n = 1:Ntex
    uuex(1,n+1) = gsin(n*dtex);
    
    for k = 2:(Nx-1)
        % a > 0, reason for Upwind
        uuex(k,n+1) = uuex(k,n) - a*lambdaex*(uuex(k,n)-uuex(k-1,n));
    end
    
    uuex(k+1,n+1) = 2*uuex(k,n+1)-uuex(k-1,n+1);
end


figure(1)
plot(xvec,uu(:,end))
hold on
plot(xvec,ulf(:,end))
hold on
plot(xvec,ulw(:,end))
hold on
plot(xvec,uuex(:,end))
legend('Upwind', 'Lax-Friedrich', 'Lax-Wendroff', 'Exact solution')
xlabel('x')
ylabel('u(x,6)')


%% 1.a) Square wave

uu = zeros(Nx,1);
ulf = uu;
ulw = uu;
uuex = uu;


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


dtex = dx;
Ntex = tend/dtex;
lambdaex = dtex/dx;

for n = 1:Ntex
    uuex(1,n+1) = gsq(n*dtex);
    
    for k = 2:(Nx-1)
        % a > 0, reason for Upwind
        uuex(k,n+1) = uuex(k,n) - a*lambdaex*(uuex(k,n)-uuex(k-1,n));
    end
    
    uuex(k+1,n+1) = 2*uuex(k,n+1)-uuex(k-1,n+1);
end



figure(2)
plot(xvec,uu(:,end))
hold on
plot(xvec,ulf(:,end))
hold on
plot(xvec,ulw(:,end))
hold on
plot(xvec,uuex(:,end))
legend('Upwind', 'Lax-Friedrich', 'Lax-Wendroff', 'Exact solution')
xlabel('x')
ylabel('u(x,6)')


%% 1.b) Sin wave
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


%% 1.b) Square wave
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


%% Part 2
clear all, close all, clc


L = 5;
v = 1;
Tcool = 50;
Thot = 200;
alpha = 0.5;

c = alpha*Tcool;
b = alpha;

a = v;
tend = 6;

Nx = 100;
dx = L/(Nx-1);
% dt <= dx
dt = 0.7*dx;
Nt = tend/dt;
lambda = dt/dx;
sigma = a*(dt/dx);

% x,t - interval
xvec = 0:dx:L;
tvec = 0:dt:tend;

% Initial condition: u(x,0) = Tcool
% Upwind: Tu; Lax-Wendroff: Tlw
Tu = Tcool*ones(Nx,1);
Tlw = Tu;


%% 2.a) 3D - plots of solution

% for d = 0.1:0.01:1
    dt = 0.7*dx;
    Nt = tend/dt;
    lambda = dt/dx;
    sigma = a*(dt/dx);
    
    % Initial condition: u(x,0) = Tcool
    Tu = Tcool*ones(Nx,1);
    Tlw = Tu;
    
    for n = 1:Nt
        Tu(1,n+1) = boundary(n*dt, Tcool, Thot);
        Tlw(1,n+1) = boundary(n*dt, Tcool, Thot);

        for k = 2:(Nx-1)
            % a > 0, reason for Upwind
            Tu(k,n+1) = Tu(k,n) - a*lambda*(Tu(k,n)-Tu(k-1,n)) - dt*(b*Tu(k,n)-c);
            Tlw(k,n+1) = Tlw(k,n) - ((sigma*(1-b*dt))/2)*(Tlw(k+1,n)-Tlw(k-1,n)) + ((sigma^2)/2)*(Tlw(k+1,n)-2*Tlw(k,n)+Tlw(k-1,n)) - dt*(1-((b*dt)/2))*(b*Tlw(k,n) - c);
            
        end

        Tu(k+1,n+1) = 2*Tu(k,n+1)-Tu(k-1,n+1);
        Tlw(k+1,n+1) = 2*Tlw(k,n+1)-Tlw(k-1,n+1);
    end

    figure1 = figure('Position', [500, 80, 600, 900])
    subplot(2,1,1)
    meshy=linspace(0,tend,1+Nt);
    meshx=linspace(0,L,Nx);
    mesh(meshy,meshx,Tu)
    title('Upwind (lambda = 0.7)')
    xlabel('t')
    ylabel('x')
    zlabel('T')
    
    subplot(2,1,2)
    mesh(meshy,meshx,Tlw)
    title('Lax-Wendroff (lambda = 0.7)')
    xlabel('t')
    ylabel('x')
    zlabel('T')
% end


%% 2.b) Different dt
close all, clc
   
for d = 0.1:0.1:1
    dt = d*dx;
    Nt = tend/dt;
    lambda = dt/dx;
    sigma = a*(dt/dx);
    meshy=linspace(0,tend,1+Nt);
    meshx=linspace(0,L,Nx);
    
    % Initial condition: u(x,0) = Tcool
    Tu = Tcool*ones(Nx,1);
    Tlw = Tu;
    
    for n = 1:Nt
        Tu(1,n+1) = boundary(n*dt, Tcool, Thot);
        Tlw(1,n+1) = boundary(n*dt, Tcool, Thot);

        for k = 2:(Nx-1)
            % a > 0, reason for Upwind
            Tu(k,n+1) = Tu(k,n) - a*lambda*(Tu(k,n)-Tu(k-1,n)) - dt*(b*Tu(k,n)-c);
            Tlw(k,n+1) = Tlw(k,n) - ((sigma*(1-b*dt))/2)*(Tlw(k+1,n)-Tlw(k-1,n)) + ((sigma^2)/2)*(Tlw(k+1,n)-2*Tlw(k,n)+Tlw(k-1,n)) - dt*(1-((b*dt)/2))*(b*Tlw(k,n) - c);

        end

        Tu(k+1,n+1) = 2*Tu(k,n+1)-Tu(k-1,n+1);
        Tlw(k+1,n+1) = 2*Tlw(k,n+1)-Tlw(k-1,n+1);
    end
   

    subplot(2,1,1)
    plot(xvec,Tu(:,end))
    hold on
    plot(xvec,Tlw(:,end))
    title('t=6')
    legend('Upwind', 'Lax-Wendroff')
    xlabel('x')
    ylabel('T')
    hold off

    subplot(2,1,2)
    plot(xvec,Tu(:,round(3/dt+1)))
    hold on
    plot(xvec,Tlw(:,round(3/dt+1)))
    title('t=3')
    legend('Upwind', 'Lax-Wendroff')
    xlabel('x')
    ylabel('T')
    drawnow
    pause(0.01)
    hold off
end

%% 2.b) Different Nx
close all, clc

for Nx = [10:20:200]
    dx = L/(Nx-1);
    % dt <= dx
    dt = 0.6*dx;
    Nt = tend/dt;
    lambda = dt/dx;
    sigma = a*(dt/dx);
    
    xvec = 0:dx:L;
    tvec = 0:dt:tend;
    
    % Initial condition: u(x,0) = Tcool
    Tu = Tcool*ones(Nx,1);
    Tlw = Tu;
    
    for n = 1:Nt
        Tu(1,n+1) = boundary(n*dt, Tcool, Thot);
        Tlw(1,n+1) = boundary(n*dt, Tcool, Thot);

        for k = 2:(Nx-1)
            % a > 0, reason for Upwind
            Tu(k,n+1) = Tu(k,n) - a*lambda*(Tu(k,n)-Tu(k-1,n)) - dt*(b*Tu(k,n)-c);
            Tlw(k,n+1) = Tlw(k,n) - ((sigma*(1-b*dt))/2)*(Tlw(k+1,n)-Tlw(k-1,n)) + ((sigma^2)/2)*(Tlw(k+1,n)-2*Tlw(k,n)+Tlw(k-1,n)) - dt*(1-((b*dt)/2))*(b*Tlw(k,n) - c);

        end

        Tu(k+1,n+1) = 2*Tu(k,n+1)-Tu(k-1,n+1);
        Tlw(k+1,n+1) = 2*Tlw(k,n+1)-Tlw(k-1,n+1);
    end
    
    
    subplot(2,1,1)
    plot(xvec,Tu(:,end))
    hold on
    plot(xvec,Tlw(:,end))
    title('t=6')
    legend('Upwind', 'Lax-Wendroff')
    xlabel('x')
    ylabel('T')
    hold off

    subplot(2,1,2)
    plot(xvec,Tu(:,round(3/dt+1)))
    hold on
    plot(xvec,Tlw(:,round(3/dt+1)))
    title('t=3')
    legend('Upwind', 'Lax-Wendroff')
    xlabel('x')
    ylabel('T')
    drawnow
    pause(0.1)
    hold off
end
