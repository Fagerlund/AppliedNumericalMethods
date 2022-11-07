function [res xk] = conjugate(A,b,eps,xk)
%GRADIENT Summary of this function goes here
%   Detailed explanation goes here

M = diag(diag(A));
M = sparse(M);
T = M-A;
T = sparse(T);

pk = xk;
rk = b;
betak = 0;

temp = 1+xk;
i = 0;

tic
while norm(rk) > eps
    pk = rk + betak*pk;
    if i > 0
        betak = (rk'*rk)/(rktemp'*rktemp);
    end
    alphak = (pk'*rk)/(pk'*A*pk);
    temp = xk;
    xk = xk + alphak*pk;
    rktemp = rk;
    rk = rk - alphak*(A*pk);
    
    
    i = i+1;
    res(i) = norm(A*xk-b)/norm(b);
end
toc

end

