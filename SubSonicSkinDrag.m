function [drag] = SubSonicSkinDrag(alt, vel, Deneb)
    [T_A, a, P_A, density_A] = atmosisa(alt, "extended","on", "action","None");
    
    mu = 1.458*10^-6*T_A^(3/2)/(T_A+110.4); % Pa*s
    
    OD = Deneb.DIAMETER; % m
    
    NC_L = Deneb.NOSECONE_LENGTH; % m
    NC_A = pi()*(OD/2)*((OD/2)+sqrt(NC_L^2+(OD/2)^2)); % m^2
    NC_SL = sqrt(NC_L^2 + (OD/2)^2); % m slant length
    [Re_NC, CD_NC, FD_NC] = skinDrag(density_A, vel, NC_SL, NC_A, mu);
    
    B_L = Deneb.BODY_LENGTH; % m body length
    B_A = pi()*OD*B_L; % m^2
    [Re_BT, CD_BT, FD_BT] = skinDrag(density_A, vel, B_L, B_A, mu);
    
    Fin_A = (Deneb.FIN.ROOT_CHORD+Deneb.FIN.TIP_CHORD)/2*Deneb.FIN.SPAN; % m^2
    Fin_L = (Deneb.FIN.ROOT_CHORD+Deneb.FIN.TIP_CHORD)/2; % m
    [Re_Fin, CD_Fin, FD_Fin] = skinDrag(density_A, vel, Fin_L, Fin_A, mu);

    drag = FD_NC + FD_BT + FD_Fin;
    if (vel == 0)
        drag = 0;
    end

    if imag(drag) ~= 0
        drag
        vel
        alt
        [Re_Fin, Re_BT, Re_NC]
    end

    function [Re,CD,FD] = skinDrag(density,V,L,A,mu)
        Re = abs(density*V*L/mu);
        
        if Re < 500000 % Turbulent
            %fprintf("Turbulent\n");
            CD = 0.031/Re^(1/7);
        elseif Re > 500000*10 % Laminar
            %fprintf("Laminar\n");
            CD = 1.328*Re^(-1/2);
        else % Mixed Flow
            %fprintf("Mixed Flow\n");
            CD = 0.031/Re^(1/7) - 1440/Re;
        end
    
        FD = CD*1/2*density*A*V^2; % N
    end
end
