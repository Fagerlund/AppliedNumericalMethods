

function out = dxODE(t, X)
    
    r1 = 0.04;
    r2 = 1e4;
    r3 = 3e7;

    x1 = X(1);
    x2 = X(2);
    x3 = X(3);

    out = [((-r1)*x1) + (r2*x2*x3); (r1*x1) - (r2*x2*x3) - (r3*(x2^2)); (r3*(x2^2))];
end

