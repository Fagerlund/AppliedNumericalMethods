function Q = Q(z)
    a = 0.2;
    b = 0.3;
    L = 1;
    Q0 = 4000;
    
    if (z >= 0 && z < a) || z > b
        Q = 0;
    else
        Q = Q0*sin((z-a)*pi/(b-a));
    end
end

