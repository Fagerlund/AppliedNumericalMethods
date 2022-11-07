clear all; close all; clc


L=1;
a=0.2;
b=0.3;
Q0=4000;
a0=100;
Tout=20;
T0=50;


h=0.05;
z=h;
n=1;
T=zeros(L/h,1);



while z<L
            
    
    
    q=Q(z);
    
    
    n=n+1;
    z=z+h
    T(n)=q;
    
end
length(T)
v=10;
A=zeros(L/h,L/h);

A= A + diag(zeros(1,L/h)+2);
A= A + diag(zeros(1,L/h-1)-1+h*v/2,1);
A= A + diag(zeros(1,L/h-1)-1-h*v/2,-1);
A=1/(h*h)*A;
Bottom=zeros(1,L/h);
Bottom(L/h)=-2*h*alfa(v);
Bottom(L/h-1)=(h*h)*A(2,1)+(h*h)*A(1,2);
Bottom=Bottom/(h*h);
A(length(A),:)=Bottom;
T(length(T))=-A(1,2)*2*h*alfa(v)*Tout;
T(1)=(1/(h*h)+v/(2*h))*T0;

aj=A(2,1);
bj=A(1,1);
cj=A(1,2);


A(L/h,L/h-1)=-bj;
A(L/h,L/h)=bj-(cj*2*h*alfa(v));


U=A\T;
U=[T0;U];
plot(linspace(0,1,L/h+1),U)


function alf=alfa(v)
a0=100;

alf=sqrt(v*v/4+a0*a0)-v/2;
end