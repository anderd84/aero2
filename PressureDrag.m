function drag_force = PressureDrag(alt, Deneb)
    [T_A, a, P_A, density_A] = atmosisa(alt, "extended","on", "action","None");

    area = Deneb.DIAMETER^2*pi/4;

    drag_force = P_A*area;
    drag_force = 0;
end