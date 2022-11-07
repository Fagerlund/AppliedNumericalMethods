function out = alph(tau)
    if tau >= 0 && tau <= 1
        out = sin(pi*tau);
    elseif tau > 1
        out = 0;
    end
end