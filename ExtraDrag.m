function drag_force = ExtraDrag(h, v, Deneb)
    [~, a, P, rho] = atmosisa(h, "extended","on", "action","None");
    mach = v/a;
    if mach < .85
        Cd = .6;
    elseif mach < 1
        Cd = .8;
    elseif mach < Deneb.LIM_TRANS
        Cd = 1.05;
    else
        Cd = .3;
    end
    drag_force = .5*Cd*rho*(Deneb.DIAMETER^2*pi/4)*v^2;
end