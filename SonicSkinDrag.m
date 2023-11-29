function drag_force = SonicSkinDrag(alt, vel, Deneb)
    R = 287.05; % J/kg*K
    k = 1.4;
    
    [T_A, a, P_A, density_A] = atmosisa(alt, "extended","on", "action","None");
    Ma_A = vel/a;
    
    OD = Deneb.DIAMETER; % m
    NC_L = Deneb.NOSECONE_LENGTH; % m


    %--------------Nosecone-------------
    NC_SL = sqrt(NC_L^2 + (OD/2)^2); % m slant length

    % oblique shock
    NC_a = asind((OD/2)/NC_SL); % deg
    [NC_WA,M_2,P_2,T_2] = obliqueShock(NC_a,Ma_A,P_A,T_A);
    
    % skin drag
    NC_A = pi()*(OD/2)*((OD/2)+sqrt(NC_L^2+(OD/2)^2)); % m^2
    
    V_2 = M_2*sqrt(k*R*T_2); % km/s
    density_2 = P_2/(R*T_2); % kg/m^3
    mu_2 = 1.458*10^-6*T_2^(3/2)/(T_2+110.4); % Pa*s
    [NC_Re,NC_CD,NC_FD] = skinDrag(density_2,V_2,NC_SL,NC_A, mu_2);

    %-------------Body-----------------
    B_L = Deneb.BODY_LENGTH; % m body length
    B_A = pi()*OD*B_L; % m^2
    
    % Prandtl Expansion
    [B_WA1,B_WA2,M_3,P_3,T_3] = prandtlExpansion(NC_a,M_2,P_2,T_2);
    
    % skin drag
    V_3 = M_3*sqrt(k*R*T_3); % km/s
    density_3 = P_3/(R*T_3); % kg/m^3
    mu_3 = 1.458*10^-6*T_3^(3/2)/(T_3+110.4); % Pa*s
    [B_Re,B_CD,B_FD] = skinDrag(density_3,V_3,B_L,B_A, mu_3);

    %---------------Fin-------------------
    Fin_T = Deneb.FIN.THICKNESS; % m
    Fin_EL = Deneb.FIN.LEADING_EDGE_PERCENT_ROOT*Deneb.FIN.ROOT_CHORD; % m leading edge length
    Fin_SL = sqrt(Fin_EL^2 + (Fin_T/2)^2); % m slant length
    
    % Oblique Shock
    Fin_a = asind((Fin_T/2)/Fin_SL); % deg
    [Fin_WA,M_4,P_4,T_4] = obliqueShock(Fin_a,M_3,P_3,T_3);
    
    % Prandtl Expansion
    [B_WA1,B_WA2,M_5,P_5,T_5] = prandtlExpansion(Fin_a,M_4,P_4,T_4);
    
    Fin_L = (Deneb.FIN.ROOT_CHORD+Deneb.FIN.TIP_CHORD)/2; % m
    Fin_A = (Deneb.FIN.ROOT_CHORD+Deneb.FIN.TIP_CHORD)/2*Deneb.FIN.SPAN; % m^2
    
    % skin drag
    V_5 = M_5*sqrt(k*R*T_5); % km/s
    density_5 = P_5/(R*T_5); % kg/m^3
    mu_5 = 1.458*10^-6*T_5^(3/2)/(T_5+110.4); % Pa*s
    [Fin_Re,Fin_CD,Fin_FD] = skinDrag(density_5,V_5,Fin_L,Fin_A, mu_5);

    drag_force = NC_FD + B_FD + Fin_FD; % N

    function [waveAngle,M2,P2,T2] = obliqueShock(d,M1,P1,T1)
        k = 1.4;
        theta = sym('theta');
        eqn = tand(d) == (2/tand(theta))*(M1^2*sind(theta)^2-1)/(M1^2*(k+cosd(2*theta))+2);
        waveAngle = double(vpasolve(eqn,theta,[0 90]));
        if sum(size(waveAngle)) < 2
            disp("EROORR!!!! no solultion for theta, are you doing subsonic?")
            waveAngle = 45;
        end
        M2 = sqrt((M1^2*sind(waveAngle)^2+2/(k-1))/((2*k)/(k-1)*M1^2*sind(waveAngle)^2-1)/sind(waveAngle-d)^2);
        P2 = P1*(2*k/(k+1)*M1^2*sind(waveAngle)^2-(k-1)/(k+1));
        T2 = T1*(7*M1^2*sind(waveAngle)^2-1)*(M1^2*sind(waveAngle)^2+5)/(36*M1^2*sind(waveAngle)^2);
    end
    
    function [waveAngle1,waveAngle2,M2,P2,T2] = prandtlExpansion(d,M1,P1,T1)
        k = 1.4;
        v1 = sqrt((k+1)/(k-1))*atand(sqrt((k-1)/(k+1)*(M1^2-1)))-atand(sqrt(M1^2-1));
        v2 = v1 + d;
    
        M2 = sym('M2');
        eqn = v2 == sqrt((k+1)/(k-1))*atand(sqrt((k-1)/(k+1)*(M2^2-1)))-atand(sqrt(M2^2-1));
        M2 = double(vpasolve(eqn,M2,[0 100]));
        if sum(size(M2)) < 2
            disp("EROORR!!!! no solultion for m2, are you doing subsonic?")
            M2 = M1;
        end
        T2 = T1*((1+(k-1)/2*M1^2)/(1+(k-1)/2*M2^2));
        P2 = P1*((1+(k-1)/2*M1^2)/(1+(k-1)/2*M2^2))^(k/(k-1));
    
        waveAngle1 = asind(1/M1);
        waveAngle2 = asind(1/M2);
    end
    
    function [Re,CD,FD] = skinDrag(density,V,L,A,mu)
        Re = density*V*L/mu;
        
        CD = 0;
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
