function [Fd, sd, pd, ed] = TotalDrag(h,v, Deneb)
    [~, a, P, ~] = atmosisa(h, "extended","on", "action","None");
    if v/a > Deneb.LIM_TRANS
        [sd, P] = SonicSkinDrag(h, v, Deneb);
    elseif v/a <=Deneb.LIM_TRANS && v/a >= 1
        sd = SubSonicSkinDrag(h, v, Deneb);
    else
        sd = SubSonicSkinDrag(h, v, Deneb);
    end
    pd = abs(PressureDrag(P, h, Deneb));
    ed = ExtraDrag(h, v, Deneb);
    Fd = sd+pd+ed;
end