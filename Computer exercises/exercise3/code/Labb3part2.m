clear all; close all ; clc
N = 10;
Tstop = 2;
Xlength = 1/2; %Given from previous task 
deltaX = (Xlength-0)/(N+1);
deltaT = 1/10*deltaX^2;


A =  diag(-2*ones(1,N)) + diag(1*ones(1,(N)-1),1) + diag(1*ones(1,(N)-1),-1);
A = 1/deltaX^2*A;
A(N+1,N)=1/deltaX^2*2;
A(N+1,N+1)=1/deltaX^2*-2;
A(N,N+1)=1/deltaX^2;
A=A*1/4
A=sparse(A);
%%
%Explicit euler, U0=0

U=zeros(1,N+1);
U=U';
AppendedArray=[];
AppendedArray=[AppendedArray,U];
tau=deltaT;
B=zeros(1,N+1);
FirstRow=[];
FirstRow=[FirstRow,0];
for i=2:Tstop/deltaT
    if tau<1;
        B(1)=1/4*sin(pi*tau)/deltaX^2;
        C=sin(pi*tau);
    else;
        B=zeros(1,N+1);
        C=0;
    end
    U=U+1/4*(deltaT*(A*U+B'));
    AppendedArray=[AppendedArray, U];

    FirstRow=[FirstRow,C];
    tau=tau+deltaT;
end


AppendedArray=[FirstRow; AppendedArray];
mesh(AppendedArray)

%% Plott
Valfor1=round(Tstop/deltaT/2);
plot(AppendedArray(:,Valfor1))

plot(AppendedArray(1,:))
hold on
plot(AppendedArray(N,:))

%% ODE solver


tspan=[0,2];
y0=zeros(1,N+1);
zeross=zeros(1,N);
Bfunc= @(t,x) [func(t,deltaX),zeross];
%opts = odeset('Jacobian',A);
tic
[t,y]= ode23(@(t,x) A*x+1/4*Bfunc(t,x)',tspan,y0);%,opts);
toc
Timesteps=length(t)
plot(t,y)











