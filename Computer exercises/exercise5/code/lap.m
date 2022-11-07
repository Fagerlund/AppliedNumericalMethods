function A = lap(n,d)
% LAP
%    A = LAP(N,D) returns the system matrix corresponding to the 
%    standard second order discretization of the Laplace operator 
%    for D dimensions and N unknowns in each coordinate direction.
%    The size of the matrix is thus N^D times N^D.

e = ones(n,1);
A1 = -spdiags([e -2*e e], -1:1, n, n);
I1 = speye(n,n);

A=A1;
I=I1;
for k=2:d
  A = kron(A,I1)+kron(I,A1);
  I = kron(I,I1);
end


