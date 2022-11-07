

function out = uvecODE(t, U)

    my = 1/82.45;
    r0 = [-my, 0]';
    r1 = [1-my, 0]';
    N = [0, 1; -1 0];
    
    u1 = U(1);
    u2 = U(2);
    u3 = U(3);
    u4 = U(4);
    
    u1u2 = (-(1-my)*(([u1; u2]-r0)./(vecnorm([u1; u2]-r0).^3))) - (my*(([u1; u2]-r1)./(vecnorm([u1; u2]-r1).^3))) + (2*N*[u3; u4]) + [u1; u2];
    out = [u3; u4; u1u2(1); u1u2(2)];
end

