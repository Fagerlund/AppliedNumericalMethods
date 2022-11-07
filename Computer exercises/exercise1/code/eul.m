

function [u1, u2] = eul(u1, u2, h, n)
    
    my = 1/82.45;
    r0 = [-my, 0]';
    r1 = [1-my, 0]';
    N = [0, 1; -1 0];
    uvec = @(u1,u2) [u2, (-(1-my)*((u1-r0)./(vecnorm(u1-r0).^3))) - (my*((u1-r1)./(vecnorm(u1-r1).^3))) + (2*N*u2) + u1];

    
    k1 = uvec(u1(:,n), u2(:,n));

    u1(:,n+1) = u1(:,n) + h*k1(:,1);
    u2(:,n+1) = u2(:,n) + h*k1(:,2);
end