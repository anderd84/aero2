function [drag] = SubSonicSkinDrag(alt, vel, Deneb)
    R = 287.05; % J/kg*K
    k = 1.4;
    
    
    [T_A, a, P_A, density_A] = atmosisa(alt, "extended","on", "action","None");
    
    CD = 0.002;
    
    OD = Deneb.DIAMETER; % m
    
    NC_L = Deneb.NOSECONE_LENGTH; % m
    NC_A = pi()*(OD/2)*((OD/2)+sqrt(NC_L^2+(OD/2)^2)); % m^2
    
    B_L = Deneb.BODY_LENGTH; % m body length
    B_A = pi()*OD*B_L; % im^2
    
    Fin_A = (Deneb.FIN.ROOT_CHORD+Deneb.FIN.TIP_CHORD)/2*Deneb.FIN.SPAN; % m^2

    drag = CD*1/2*density_A*((NC_A+B_A+Fin_A))*vel^2; % N
end
