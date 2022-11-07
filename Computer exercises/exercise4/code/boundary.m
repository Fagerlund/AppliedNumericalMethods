function T0 = boundary(t, Tcool, Thot)
    if t <= 0.125
        T0 = Tcool + (Thot - Tcool)*sin(4*pi*t);
    elseif t > 0.125 && t <= 1
        T0 = Thot;
    elseif t > 1
        T0 = Thot + Tcool*sin(5*pi*(t-1));
    end
end

