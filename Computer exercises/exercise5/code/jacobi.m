function [res xk] = jacobi(A,b,eps,xk)
%JACOBI Summary of this function goes here
%   Detailed explanation goes here

M = diag(diag(A));
M = sparse(M);
T = M-A;
T = sparse(T);

% to start the while loop
temp = 1+xk;
i = 0;

tic
while norm(xk-temp) > eps
    temp = xk;
    xk = M\(T*xk+b);
    
    r = xk - temp;
    i = i+1;
    res(i) = norm(A*xk-b)/norm(b);
end
toc

end
