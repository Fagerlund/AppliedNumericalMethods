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


%% a) 3D - plots of solution

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
    meshx=linspace(0,L,Nx)
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


%% b) Different dt
close all, clc
   
for d = 0.1:0.01:1
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

%% Different Nx
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