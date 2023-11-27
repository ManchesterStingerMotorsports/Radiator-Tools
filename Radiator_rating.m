% Radiator Rating Tool
% Tom Wilkinson tomwilkinson7@me.com
%% THIS HAS NOT BEEN EXPERIMENTALY VERIFIED %%

clearvars; format shortG

% Modify these

% Radiator Specs
radiator_L = 0.605; % [m] 
radiator_W = 0.058; % [m] 
radiator_H = 0.145; % [m] 
tubes_N = 17; 
tube_H = 0.0015; % [m] 
fin_spacing = 25.4/18000; % [m]   % if the fins are /\/\/\/\
wall_thickness = 1e-4; % [m]   /\ = 1 fin spacing


% Flow Rates
water_flowrate = 5.5; % [L/min]
air_flowrate = .1; % [m^3/s]
water_temp = 32; % [degC]
air_temp = 22; % [degC]

%% DO NOT EDIT ANYTHING BELOW HERE %% 

% Derived specs
tube_H_internal = tube_H - 2*wall_thickness; % [m] 
tube_W_internal = radiator_W - 2*wall_thickness; % [m] 
gap_H = (radiator_H - tubes_N * tube_H) / (tubes_N); % [m] Fin height 
fins_N = tubes_N * radiator_L / fin_spacing; 
flow_speed = (water_flowrate/60000)/(tube_W_internal*tube_H_internal*tubes_N); 
air_speed = ((air_flowrate) / (radiator_L * radiator_H)) ; 

A_c = (2*fin_spacing + 4*((fin_spacing/2)^2 + gap_H^2)^.5)*radiator_W  * fins_N; 
A_h = 2 * (tube_H_internal + tube_W_internal) * fin_spacing  * fins_N; 

% Uses data from Compact Heat Exchangers Third Edition, Kays and London
% (and a lot of linear interpolation)

% Air side - Surface 14.77
c_pc = interp1([0, 10, 20, 30, 40], [1004, 1004, 1004, 1005, 1005], air_temp, "linear", "extrap");
u_c = interp1([0, 10, 20, 30, 40], [17.2e-6, 17.69e-6, 18.17e-6, 18.64e-6, 19.11e-6], air_temp, "linear", "extrap"); 
Pr_c = interp1([0, 10, 20, 30, 40], [.717, .714, .712, .710, .709], air_temp, "linear", "extrap");
r_hc = 6.475e-4; 
rho_c = interp1([0, 10, 20, 30, 40], [1.293, 1.2473, 1.2047, 1.165, 1.1277], air_temp, "linear", "extrap");
Re_c = (rho_c * air_speed * 4 * r_hc)/u_c;
StPr_c = interp1([10000, 8000, 6000, 5000, 4000, 3000, 25000, 2000, 1500, 1200, 1000, 800, 600, 500], [.00310, .00326, .00352, .00367, .00389, .00417, .00435, .00456, .00495, .00538, .00585, .00663, .00791, .00898], Re_c, "linear", "extrap");

h_c = (c_pc * u_c * StPr_c * Re_c)/((Pr_c^(2/3)) * 4 * r_hc); 


% Water side - Surface FT-1
c_ph = interp1([0.01, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 140], [4217, 4193, 4182, 4179, 4179, 4181, 4185, 4190, 4197, 4205, 4216, 4285], water_temp, "linear", "extrap"); 
u_h = interp1([0.01, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 140], [17.91e-4, 13.08e-4, 10.03e-4, 79.77e-5, 65.31e-5, 54.71e-5, 46.68e-5, 40.44e-5, 35.49e-5, 31.50e-5, 28.22e-5, 19.61e-5], water_temp, "linear", "extrap"); 
Pr_h = interp1([0.01, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 140], [13.44, 9.42, 6.99, 5.42, 4.34, 3.57, 3.00, 2.57, 2.23, 1.97, 1.76, 1.23], water_temp, "linear", "extrap");
r_hh = (tube_W_internal*tube_H_internal)/(2*(tube_W_internal+tube_H_internal));
rho_h = interp1([0.01, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 140], [999.8, 999.8, 998.2, 995.6, 992.2, 998.0, 983.2, 977.7, 971.8, 965.3, 958.3, 926.1], water_temp, "linear", "extrap");
Re_h = (rho_h * flow_speed * 4 * r_hh)/(u_h);
StPr_h = interp1([10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1200, 1000, 800, 600, 500], [.00292, .00298, .00305, .00310, .00317, .00291, .00272, .00310, .00381, .00453, .00522, .00621, .00780, .00900], water_temp, "linear", "extrap");

h_h = (c_ph * u_h * StPr_h * Re_h)/((Pr_h^(2/3)) * 4 * r_hh); 

% Fin efficiency
m = ((2 * h_c)/(240 * wall_thickness))^.5;
l = ((fin_spacing/2)^2 + gap_H^2)^.5 / 2;
n_f = tanh(m*l)/(m*l);
n_c = 1 - (4*((fin_spacing/2)^2 + gap_H^2)^.5/(tube_H * fin_spacing + 2*fin_spacing + 4*((fin_spacing/2)^2 + gap_H^2)^.5)) * (1-n_f);
% Overall heat transfer coefficients
U_h = 1/((1/h_h) + (1/((A_c/A_h)*h_c*n_c))); 
U_c = 1/(1/(h_c*n_c) + 1/((A_h/A_c)*h_h)); 

% E-NTU Method
UA = U_c * A_c; 
c_max = max((water_flowrate/60000) * rho_h * c_ph, air_flowrate * rho_c * c_pc); 
c_min = min((water_flowrate/60000) * rho_h * c_ph, air_flowrate * rho_c * c_pc); 

c_r = c_min / c_max;
NTU = UA / c_min;

effectiveness = interp2([0, 0.25, 0.5, 0.75, 1], ... % Crossflow heat exchanger, both fluids unmixed
    [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7], ...
    [0,0,0,0,0; ...
    .221, .215, .209, .204, .199; ...
    .393, .375, .358, .341, .326; ...
    .528, .495, .466, .439, .413; ...
    .632, .588, .547, .510, .476; ...
    .714, .660, .610, .565, .523; ...
    .777, .716, .660, .608, .560; ...
    .826, .761, .700, .642, .590; ...
    .865, .797, .732, .671, .614; ...
    .918, .851, .783, .716, .652; ...
    .950, .888, .819, .749, .681; ...
    .970, .915, .848, .776, .704; ...
    .982, .934, .869, .797, .722; ...
    .989, .948, .887, .814, .737; ...
    .993, .959, .901, .829, .751; ...
    .997, .974, .924, .853, .772; ...
    .999, .983, .940, .871, .789], c_r, NTU, 'spline');

disp(['U_c: ' num2str(U_c)]) % U_c should be between 25 and 50 (ish)
disp(['Heat dissipation: ' num2str(round(effectiveness * c_min * (water_temp - air_temp))) 'W'])