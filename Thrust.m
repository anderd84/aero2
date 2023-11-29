function [thrust, mdot] = Thrust(t, P_atm, Deneb)
    % R_universal = 1545.349; % lbf*ft/(lbmol*R)
    % Molar_Mass = 22.54; % (lbm/lbmol) (From RPA)
    % R_gas = R_universal/Molar_Mass; % lbf*ft/(lbm*R)
    % T_Chamber = 6229.0292; %(R)
    % g_c = 32.174; %(lbm ft/ lbf s^2)
    % gamma = 1.1787;
    % A_Exit = pi .* (1/4) .*(Deneb.ENGINE.EXIT_AREA).^2; %(in^2)
    % 
    % d_throat = linspace(d_throat_nominal, d_throat_nominal+BurnTime * Throat_ablation/1000, index); %(in)
    % A_throat = pi .* (1/4) .*(d_throat).^2; %(in^2)
    % P_Chamber = mdot_calc .* (sqrt(g_c * R_gas * T_Chamber))./ A_throat ./ sqrt(gamma) ./((gamma +1)  / 2)^(-(gamma +1)/(2* (gamma -1)));

    thrust = 2000 * 4.44822; % N
    mdot = 7.5 / 2.20462;

    % function [Exit_Mach, Exit_Pressure, Exit_Temperature] = Nozzle(Total_Temperature, Total_Pressure, Throat_Area, Exit_Area, Specific_Heat_Ratio)
    %     syms Ma;
    %     Ma_exit = vpasolve(Exit_Area/Throat_Area == (1/Ma) *((1+((Specific_Heat_Ratio-1)/2)*(Ma^2))/(1+(Specific_Heat_Ratio-1)/2))^((Specific_Heat_Ratio+1)/(2*(Specific_Heat_Ratio-1))), Ma, 3);
    %     if(length(Ma_exit) <2)
    %         Exit_Mach = double(Ma_exit(1));
    %     else 
    %         Exit_Mach = double(Ma_exit(2));
    %     end
    %     Exit_Pressure = double(Total_Pressure/(1+((Specific_Heat_Ratio-1)/2)*Exit_Mach^2)^(Specific_Heat_Ratio/(Specific_Heat_Ratio-1)));
    %     Exit_Temperature = double(Total_Temperature/(1+((Specific_Heat_Ratio-1)/2)*Exit_Mach^2));
    % end

end

