%function q = Rate(t, s)
function q = Rate(car_speed, inlet_temp)
% Radiator Specs
radiator_L = 0.304; % [m]
radiator_W = 0.0399; % [m]
radiator_H = 0.208; % [m]
tubes_N = 20;
tube_H = 0.0015; % [m]
fin_spacing = 0.002; % [m]
wall_thickness = 1e-4; % [m]


% Flow Rates
water_flowrate = 20; % [L/min]
V_car = car_speed; % [m/s] car velocity
water_temp = inlet_temp; % [degC]
air_temp = 20; % [degC]

%% DO NOT EDIT ANYTHING BELOW HERE %% 

% Abandon hope ye who enter here %

% Derived specs
tube_H_internal = tube_H - 2*wall_thickness; % [m] 
tube_W_internal = radiator_W - 2*wall_thickness; % [m] 
gap_H = (radiator_H - tubes_N * tube_H) / (tubes_N); % [m] Gap between adjacent fins
fins_N = tubes_N * radiator_L / fin_spacing;

flow_speed = (water_flowrate/60000)/(tube_W_internal*tube_H_internal*tubes_N);
rho_c = interp1([0, 10, 20, 30, 40], [1.293, 1.2473, 1.2047, 1.165, 1.1277], air_temp, "linear", "extrap");
air_flowrate = Vdot(radiator_H*radiator_L, rho_c, V_car); % [m^3/s]
air_speed = ((air_flowrate) / (radiator_L * radiator_H)); 

% The surface area associated with one fin spacing
A_c = (4*((fin_spacing/2)^2 + gap_H^2)^.5)*radiator_W * fins_N + (tube_H + radiator_W * radiator_L * tubes_N * 2) ; % Air side
A_h = 2 * (tube_H_internal + tube_W_internal) * radiator_L * tubes_N; % Water side

% Uses data from Compact Heat Exchangers Third Edition, Kays and London
% (and a lot of linear interpolation)

% Air side - Surface 14.77
c_pc = interp1([0, 10, 20, 30, 40], [1004, 1004, 1004, 1005, 1005], air_temp, "linear", "extrap");
u_c = interp1([0, 10, 20, 30, 40], [17.2e-6, 17.69e-6, 18.17e-6, 18.64e-6, 19.11e-6], air_temp, "linear", "extrap"); 
Pr_c = interp1([0, 10, 20, 30, 40], [.717, .714, .712, .710, .709], air_temp, "linear", "extrap");
r_hc = 6.475e-4; % TODO would like to make this dynamic, but struggling to get the equation right
Re_c = (rho_c * air_speed * 4 * r_hc)/u_c;
StPr_c = interp1([10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1200, 1000, 800, 600, 500], [.00310, .00326, .00352, .00367, .00389, .00417, .00435, .00456, .00495, .00538, .00585, .00663, .00791, .00898], Re_c, "linear", "extrap");

h_c = (c_pc * u_c * StPr_c * Re_c)/((Pr_c^(2/3)) * 4 * r_hc); % Equation [1-1]


% Water side - Surface FT-1
c_ph = interp1([0.01, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 140], [4217, 4193, 4182, 4179, 4179, 4181, 4185, 4190, 4197, 4205, 4216, 4285], water_temp, "linear", "extrap"); 
u_h = interp1([0.01, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 140], [17.91e-4, 13.08e-4, 10.03e-4, 79.77e-5, 65.31e-5, 54.71e-5, 46.68e-5, 40.44e-5, 35.49e-5, 31.50e-5, 28.22e-5, 19.61e-5], water_temp, "linear", "extrap"); 
Pr_h = interp1([0.01, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 140], [13.44, 9.42, 6.99, 5.42, 4.34, 3.57, 3.00, 2.57, 2.23, 1.97, 1.76, 1.23], water_temp, "linear", "extrap");
r_hh = (tube_W_internal*tube_H_internal)/(2*(tube_W_internal+tube_H_internal));
rho_h = interp1([0.01, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 140], [999.8, 999.8, 998.2, 995.6, 992.2, 998.0, 983.2, 977.7, 971.8, 965.3, 958.3, 926.1], water_temp, "linear", "extrap");
Re_h = (rho_h * flow_speed * 4 * r_hh)/(u_h);
StPr_h = interp1([10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1200, 1000, 800, 600, 500], [.00292, .00298, .00305, .00310, .00317, .00291, .00272, .00310, .00381, .00453, .00522, .00621, .00780, .00900], water_temp, "linear", "extrap");

h_h = (c_ph * u_h * StPr_h * Re_h)/((Pr_h^(2/3)) * 4 * r_hh); % Equation [1-1]

% Fin efficiency
m = ((2 * h_c)/(240 * wall_thickness))^.5; % Equation [2-5]
l = ((fin_spacing/2)^2 + gap_H^2)^.5 / 2; % Fin length
n_f = tanh(m*l)/(m*l); % Equation [2-4]
n_c = 1 - (4*((fin_spacing/2)^2 + gap_H^2)^.5/(tube_H * fin_spacing + 2*fin_spacing + 4*((fin_spacing/2)^2 + gap_H^2)^.5)) * (1-n_f); % Equation [2-3]
% Overall heat transfer coefficients, Equation [2-2] 
U_h = 1/((1/h_h) + (1/((A_c/A_h)*h_c*n_c))); % Water side
U_c = 1/(1/(h_c*n_c) + 1/((A_h/A_c)*h_h)); % Air side

% E-NTU Method
UA = U_c * A_c; % U_c * A_c == U_h * A_h
c_max = max((water_flowrate/60000) * rho_h * c_ph, air_flowrate * rho_c * c_pc); 
c_min = min((water_flowrate/60000) * rho_h * c_ph, air_flowrate * rho_c * c_pc); 

c_r = c_min / c_max;
NTU = UA / c_min;
effectiveness = 1 - exp(NTU^0.22 * (exp(-c_r*NTU^0.78)-1) /c_r); % https://www.semanticscholar.org/reader/0b0893aa2fad9bbd160f2fa2e5d1aa5c92eb6b7b


q = effectiveness * c_min * (water_temp - air_temp);
end

function Q_F = Vdot(radiator_area, rho, V_car)
    % air flow rate calculation THIS PAPER IS FUCKING MAGIC
    % https://ideaexchange.uakron.edu/cgi/viewcontent.cgi?article=1143&context=honors_research_projects  
    % SPAL 190mm diameter 24V brushed fans
    fan_D = 0.190;                  % Fan diameter [m]
    n = 2;                          % number of fans
    c1  = 350/(730/3600);           % Fan coeff (slope of linear regression)
    kr  = 15.79;                    % Radiator pressure coeff (from paper, a bit of a guess)
    a1  = radiator_area;            % Radiator frontal area
    a4  = n * pi * fan_D^2 /4;      % Fan frontal area
    c0  = 350;                      % Fan coeff (intercept of linear regression)
    v   = V_car;                    % Velocity of car

    a = rho * 0.5 * ((kr/(a1^2)) + (1/(a4^2)));
    b = c1/n;
    c = -(c0 + (0.5*rho*v^2));
    Q_F = (-b + sqrt(b^2 - 4*a*c))/(2*a);
end