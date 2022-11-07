

%% 1.a Is the convergence faster or slower with larger n?
clc, close all

d = 2;
eps = 1e-4;

for n = [10 20 40 80]
    N = n^d;
    xk = zeros(N,1);
    b = rand(N,1);
    A = lap(n,d);
    A = sparse(A);


    L=["Time Jacobi, n=",n,'d=',d];
    disp(L)
    [resJ xendJ] = jacobi(A,b,eps,xk);

    L=["Time Conjugate, n=",n,'d=',d];
    disp(L)
    [resC xendC] = conjugate(A,b,eps,xk);


    semilogy(0:length(resJ)-1,resJ)
    hold on
    semilogy(0:length(resC)-1,resC)
    legend('Jacobi','Conjugate')
    title("Relative residual")
    hold off

end



%% 1.a Which method converges faster?
clear all, close all, clc
n = 40;
d = 2;
eps = 1e-4;

N = n^d;
xk = zeros(N,1);
b = rand(N,1);
A = lap(n,d);
A = sparse(A);


L=["Time Jacobi, n=",n,'d=',d];
disp(L)
[resJ xendJ] = jacobi(A,b,eps,xk);

L=["Time Conjugate, n=",n,'d=',d];
disp(L)
[resC xendC] = conjugate(A,b,eps,xk);



semilogy(0:length(resJ)-1,resJ)
hold on
semilogy(0:length(resC)-1,resC)
legend('Jacobi','Conjugate')
title("Relative residual")
hold off


%% Change number of iterations
clear all, close all, clc
% number of iterations
numb = 2000;
iter = 1:numb;

% Jacobi

n = 20;
d = 3;

N = n^d;
xk = zeros(N,1);
b = rand(N,1);
A = lap(n,d);
A = sparse(A);

M = diag(diag(A));
M = sparse(M);
T = M-A;
T = sparse(T);

lim = 50;

resJlog = [];
iterlog = [];
k = 2;
tic
for i = iter
    temp = xk;
    xk = M\(T*xk+b);
    
    r = xk - temp;
    resJ(i) = norm(A*xk-b)/norm(b);
    
        
    if i > lim
        if resJ(i) < 0.1*resJ(lim)
            break
        end
    end
end
toc



% Conjugate

N = n^d;
xk = zeros(N,1);
b = rand(N,1);
A = lap(n,d);
A = sparse(A);

M = diag(diag(A));
M = sparse(M);
T = M-A;
T = sparse(T);

pk = xk;
rk = b;
betak = 0;

resClog = [];
iterlog = [];
k = 2;
tic
for i = iter
    pk = rk + betak*pk;
    if i > 1
        betak = (rk'*rk)/(rktemp'*rktemp);
    end
    alphak = (pk'*rk)/(pk'*A*pk);
    temp = xk;
    xk = xk + alphak*pk;
    rktemp = rk;
    rk = rk - alphak*(A*pk);
    
    
    resC(i) = norm(A*xk-b)/norm(b);
    
    if i > lim
        if resC(i) < 0.1*resC(lim)
            break 
        end
    end
end
toc

semilogy(0:length(resJ)-1,resJ)
hold on
semilogy(0:length(resC)-1,resC)
legend('Jacobi','Conjugate')
title("Relative residual")
hold off

length(resJ)-lim
length(resC)-lim


%% 1.a Is the convergence faster or slower with larger d, for a fixed matrix size N?

clear all, clc, close all

eps = 1e-4;


for d = [1 2 3 4]
    Nfix = 1000;
    n = round(nthroot(Nfix,d));
    N = n^d;
    xk = zeros(N,1);
    b = rand(N,1);
    A = lap(n,d);
    A = sparse(A);

    L=["Time Jacobi, n=",n,'d=',d];
    disp(L)
    [resJ xendJ] = jacobi(A,b,eps,xk);

    L=["Time Conjugate, n=",n,'d=',d];
    disp(L)
    [resC xendC] = conjugate(A,b,eps,xk);


    semilogy(0:length(resJ)-1,resJ)
    hold on
    semilogy(0:length(resC)-1,resC)
    legend('Jacobi','Conjugate')
    title("Relative residual")
    hold off

end


%% 1.b Varying d
clear all, close all, clc

n = 10;
eps = 1e-10;

for d = 1:5
    N = n^d;
    xk = zeros(N,1);
    b = rand(N,1);
    A = lap(n,d);
    A = sparse(A);
    
    L=["Time Conjugate, n=",n,'d=',d];
    disp(L)
    [resC xendC] = conjugate(A,b,eps,xk);
    
    L = ["Time for backslash"];
    disp(L)
    tic
    xendB = A\b;
    toc
    norm(xendB-xendC)
end


%% 1.b Varying n
clear all, close all, clc

d = 3;
eps = 1e-10;

for n = 1:5:26
    N = n^d;
    xk = zeros(N,1);
    b = rand(N,1);
    A = lap(n,d);
    A = sparse(A);
    
    L=["Time Conjugate, n=",n,'d=',d];
    disp(L)
    [resC xendC] = conjugate(A,b,eps,xk);
    
    L = ["Time for backslash"];
    disp(L)
    tic
    xendB = A\b;
    toc
    norm(xendB-xendC)
end


%% Part 2

clear all;close all;clc
load cooling_flange.mat
spy(A)
b=rand(length(A),1);

%% a
clc
tol=1e-4;
tic
[X,FLAG,RELRES,ITER,RESVEC]=pcg(A,b,tol,10000);
FLAG
L=["Computational time for pcg"];
disp(L)
toc
L=["Number of iterations used", ITER];
disp(L)
L=["Final size of the relative residual", RELRES];
disp(L)
l=linspace(0,ITER,ITER+1);
semilogy(l,RESVEC)
tic
l=A\b;
L=["Computational time for Backslash"];
disp(L)
toc
xlabel("# of iterations")
ylabel("residual")





%% b
set(gca, 'YScale', 'log')
clc;
tol=1e-4;
M=diag(A);
M=diag(M);
tic
L=["Computational time for pcg Diag"];

[X,FLAG,RELRES,ITER,RESVEC]=pcg(A,b,tol,10000,M);

disp(L)
toc
L=["Number of iterations used", ITER];
disp(L)
L=["Final size of the relative residual", RELRES];
disp(L)
l=linspace(0,ITER,ITER+1);
hold on
semilogy(l,RESVEC)


l=ichol(A);
L=["Computational time for pcg Cholesky"];

lt=l';
tic

[X,FLAG,RELRES,ITER,RESVEC]=pcg(A,b,tol,10000,l,lt);

disp(L)
toc
L=["Number of iterations used", ITER];
disp(L)
L=["Final size of the relative residual", RELRES];
disp(L)

l=linspace(0,ITER,ITER+1);

semilogy(l,RESVEC)
ylabel("residual")
xlabel("# of iterations" )



tic
l=A\b;
L=["Computational time for Backslash"];
disp(L)
toc

%% c
clear all; close all; clc
load convdiff.mat
spy(A)
b=rand(length(A),1);
tol=1e-4;
L=["Computational time for GMRES"];
[X,FLAG,RELRES,ITER,RESVEC]=pcg(A,b,tol,10000); %Just verify by FLAG for high iter?
FLAG


[L1,U]=ilu(A);
tic
[X,FLAG,RELRES,ITER,RESVEC1]=gmres(A,b,[] ,tol,20000,L1,U);
disp(L)
toc
L=["Number of iterations used", ITER(2)];
disp(L)
L=["Final size of the relative residual", RELRES];
disp(L)

tic
l=A\b;
L=["Computational time for Backslash"];
disp(L)
toc


% Not enough memory, this is the result we obtained from a better computer:
% FLAG =
% 
%      4
% 
% Computational time for GMRES
% Elapsed time is 0.618964 seconds.
%     "Number of iterations used"    "70"
% 
%     "Final size of the relative residual"    "8.6606e-05"
% 
% Computational time for Backslash
% Elapsed time is 3.086448 seconds.

