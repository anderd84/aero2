function drag_force = PressureDrag(Ptip, alt, Deneb)
    area = Deneb.DIAMETER^2*pi/4;
    [~, ~, Pbase, ~] = atmosisa(alt-Deneb.BODY_LENGTH, "extended", "on", "action", "None");
    drag_force = (Ptip-Pbase)*area;
end