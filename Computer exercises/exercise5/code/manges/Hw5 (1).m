clear all; close all; clc


%% Plot relative
%Constants picked arbitrarily
n=25;
d=3;
error=0.0001;

figure(1)
resj=Jacobi(n,d,error);
xJ=linspace(0,1,length(resj));
semilogy(xJ,resj)
hold on

[resC,b,sol]=Conjugate(n,d,error);
xC=linspace(0,1,length(resC));
semilogy(xC,resC)
hold on
legend('Jacobi','Conjugate')
title("Relative residual")
hold off



%% Plot absolute
n=25;
d=3;
error=0.0001;

figure(1)
resj=Jacobi(n,d,error);
xJ=linspace(0,length(resj),length(resj));
semilogy(xJ,resj)
hold on

[resC,b,sol]=Conjugate(n,d,error);
xC=linspace(0,length(resC),length(resC));
semilogy(xC,resC)
hold on
legend('Jacobi','Conjugate')
title("Iteration curve")
hold off


%% Comparing various speeds with n
clear all;close all;clc
%Also used for determining which method converges faster
%Jacobi
%constants are set arbitrarily
d=3;
epsilon=0.0001;
ItersJ=[];
spacingN=linspace(5,29,7);
for n=5:4:29
    J=Jacobi(n,d,epsilon);
    ItersJ=[ItersJ length(J)];
    
    
end

% Conjugate
ItersC=[];
for n=5:4:29
    [C,b,sol]=Conjugate(n,d,epsilon);
    ItersC=[ItersC length(C)];
    
    
end

plot(spacingN,ItersC)
hold on 
plot(spacingN,ItersJ)



%% Comparing various speeds with d, with fixed N 
clear all;close all;clc
%Also used for determining which method converges faster
%Jacobi
%constants are set arbitrarily
N=10000;
epsilon=0.0001;
ItersJ=[];
spacingd=linspace(3,7,7-3+1);
for d=3:7
    n=round(nthroot(N,d));
    J=Jacobi(n,d,epsilon);
    ItersJ=[ItersJ length(J)];
    
    
end

% Conjugate
ItersC=[];
for d=3:7
    n=round(nthroot(N,d));
    [C,b,sol]=Conjugate(n,d,epsilon);
    ItersC=[ItersC length(C)];
    
    
end

plot(spacingd,ItersC)
hold on 
plot(spacingd,ItersJ)


%% b part
clear all; close all; clc

%Varying d
n=5;
error=1e-10;
for d=1:6
    [a,b,sol]=Conjugate(n,d,error);
    A=lap(n,d);
    L=["Time for backslash"];
    disp(L)
    tic
    sol1=A\b;
    toc
    norm(sol1-sol)
end
    


%% b part
clear all; close all; clc

%Varying n
d=3;
error=1e-10;
for n=1:5:26
    [C,b,sol]=Conjugate(n,d,error);
    A=lap(n,d);
    L=["Time for backslash"];
    disp(L)
    tic
    sol1=A\b;
    toc
    norm(sol1-sol)
end
