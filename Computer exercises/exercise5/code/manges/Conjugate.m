function [resC,b,sol] = Conjugate(n,d,epsilon)

A=lap(n,d);
N=n^d;
M=tril(A);
T=M-A;
b=rand(N,1);

%Iterating
bool=true; 
%beta = norm(inv(M)*(A-M))
x0=zeros(N,1);
xk=x0;
iterC=0;
tic
while bool==true
    
    LHS=T*xk+b;
    temp=xk;
    xk=M\(LHS);
    
    residual=temp-xk;
    resC(iterC+1)=norm(A*xk-b)/norm(b);
    iterC=iterC+1;
    if norm(residual)<epsilon
        break
    end
    sol=xk;
end
L=["Time elapsed for Conjugate, with n=",n,'d=',d];
disp(L) 
toc

end

