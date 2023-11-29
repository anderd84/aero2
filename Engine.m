function [Thrust_Newtons_tot, mdot_kg_sec, Thrust_Newtons_Velocity, Thrust_Newtons_Pressure] = Engine(timestep, P, Deneb) %all in metric

    R_universal = 1545.349; % lbf*ft/(lbmol*R)
    Molar_Mass = 22.54; %(lbm/lbmol) (From RPA)
    R_gas = R_universal/Molar_Mass; % lbf*ft/(lbm*R);
    BurnTime = 20; %(s)
    P_Atmospheric = P * 0.000145038; %psi;
    Throat_ablation = Deneb.ENGINE.THROAT_REGRESSION_RATE; %(thou/s)
    d_throat_nominal = Deneb.ENGINE.THROAT_DIAMETER_0;
    mdot_calc_lbm = Deneb.ENGINE.MDOT_0;
    mdot_calc = mdot_calc_lbm ./ 32.174; %lbm/s;
    
    
    T_Chamber = 6229.0292; %(R)
    g_c = 32.174; %(lbm ft/ lbf s^2)
    gamma = 1.1787;
    A_Exit = Deneb.ENGINE.EXIT_AREA; %(in^2)
    
    
    d_throat = d_throat_nominal + Throat_ablation *timestep;  %(in);
    A_throat = pi .* (1/4) .*(d_throat).^2; %(in^2)
    P_Chamber = mdot_calc .* (sqrt(g_c * R_gas * T_Chamber))./ A_throat ./ sqrt(gamma) ./((gamma +1)  / 2)^(-(gamma +1)/(2* (gamma -1)));
    R_gas_RPA = R_gas/778.169; % BTU/(lbm*R);
    [Exit_Mach, Exit_Pressure, Exit_Temperature] = Nozzle(T_Chamber, P_Chamber, mdot_calc_lbm,A_throat,A_Exit, gamma);
    Sonic_Velocity = sqrt(gamma * R_gas * Exit_Temperature * g_c); %ft/s;
    Exit_Velocity = Exit_Mach * Sonic_Velocity ;%ft/s;
    mdot_kg_sec = mdot_calc_lbm * 0.453592;
    Thrust_lbf_Velocity = mdot_calc_lbm /g_c * Exit_Velocity;
    Thrust_lbf_Pressure = ((Exit_Pressure) - P_Atmospheric) *A_Exit;
    Thrust_lbf_tot = Thrust_lbf_Velocity + Thrust_lbf_Pressure;
    Thrust_Newtons_Velocity = Thrust_lbf_Velocity * 4.44822;
    Thrust_Newtons_Pressure = Thrust_lbf_Pressure * 4.44822;
    Thrust_Newtons_tot = Thrust_lbf_tot * 4.44822;
    
    function [Exit_Mach, Exit_Pressure, Exit_Temperature] = Nozzle(Total_Temperature, Total_Pressure, mdot_tot,Throat_Area,Exit_Area, Specific_Heat_Ratio)
        Ma = sym('Ma');
        Ma_exit = vpasolve(Exit_Area/Throat_Area == (1/Ma) *((1+((Specific_Heat_Ratio-1)/2)*(Ma^2))/(1+(Specific_Heat_Ratio-1)/2))^((Specific_Heat_Ratio+1)/(2*(Specific_Heat_Ratio-1))), Ma, 3);
        if(length(Ma_exit) <2)
            Exit_Mach = double(Ma_exit(1));
        else 
            Exit_Mach = double(Ma_exit(2));
        end
        Exit_Pressure = double(Total_Pressure/(1+((Specific_Heat_Ratio-1)/2)*Exit_Mach^2)^(Specific_Heat_Ratio/(Specific_Heat_Ratio-1)));
        Exit_Temperature = double(Total_Temperature/(1+((Specific_Heat_Ratio-1)/2)*Exit_Mach^2));
        
    end

end
