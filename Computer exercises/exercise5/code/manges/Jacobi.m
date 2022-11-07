function resJ = Jacobi(n,d,epsilon)


A=lap(n,d);
N=n^d;
M=diag(A);
M=diag(M);
T=M-A;

b=rand(N,1);

%Iterating
bool=true; 
%beta = norm(inv(M)*(A-M))
x0=zeros(N,1);
xk=x0;
iterJ=0;
%res=zeros(3*N,1);
tic
while bool==true
    
    LHS=T*xk+b;
    temp=xk;
    xk=M\(LHS);
    
    
    residual=norm(temp-xk);
    resJ(iterJ+1)=norm(A*xk-b)/norm(b);
    iterJ=iterJ+1;
    if norm(residual)<epsilon
        break
    end
    
end
L=["Time elapsed for Jacobi, with n=",n,'d=',d];
disp(L) 
toc

end

