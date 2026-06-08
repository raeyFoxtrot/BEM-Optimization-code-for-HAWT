clc; clear; close all;

% Inputs
U = 18;                 % freestream velocity (m/s)
R = 50;                 % blade radius (m)
B = 3;                  % no. of blades
rho = 1.225;            % air density
mu = 1.8e-5;            % Dynamic viscosity

theta_root = deg2rad(0);
theta_tip  = deg2rad(5);

tol = 1e-6;
max_iterations = 500;

% Lambda sweep
lambda_range = 0:0.2:10;
Cp_list = zeros(size(lambda_range));

best_lambda = 0;
best_Cp = -inf;
best_results = struct();

% DISCRETIZATION
N = R;
r_hub = 0.2 * R;
r_vec = linspace(r_hub, R, N);
dr = r_vec(2) - r_vec(1);

% Geometry
chord = 2*(1 - r_vec/R) + 1;
theta = theta_root + (theta_tip - theta_root) * (r_vec / R);

% RE Data
Re_list = [50e3, 100e3, 200e3, 500e3, 1000e3];

polar{1} = readtable('NACA0018_50k.csv');
polar{2} = readtable('NACA0018_100k.csv');
polar{3} = readtable('NACA0018_200k.csv');
polar{4} = readtable('NACA0018_500k.csv');
polar{5} = readtable('NACA0018_1000k.csv');

% Sweep loop
for L = 1:length(lambda_range)

    lambda = lambda_range(L);
    Omega = lambda * U / R;

    a = zeros(1,N);
    ap = zeros(1,N);
    phi_arr = zeros(1,N);
    aoa_arr = zeros(1,N);
    dT = zeros(1,N);

    % BEM Loop
    for j = 1:N
        
        r = r_vec(j);
        sigma = (B * chord(j)) / (2*pi*r);
        
        a_j = 0.3;
        ap_j = 0.01;
        
        for iter = 1:max_iterations
            
            phi = atan( (U*(1 - a_j)) / (Omega*r*(1 + ap_j)) );
            
            alpha = phi - theta(j);
            alpha_deg = rad2deg(alpha);
            
            Vrel = sqrt((U*(1-a_j))^2 + (Omega*r*(1+ap_j))^2);
            Re = rho * Vrel * chord(j) / mu;
            
            idx = find(Re_list <= Re, 1, 'last');
            if isempty(idx)
                idx = 1;
            elseif idx == length(Re_list)
                idx = length(Re_list)-1;
            end
            
            Re_low = Re_list(idx);
            Re_high = Re_list(idx+1);
            
            w = (Re - Re_low) / (Re_high - Re_low);
            w = max(0, min(1, w));
            
            alpha_low = polar{idx}.Alpha;
            Cl_low_data = polar{idx}.Cl;
            Cd_low_data = polar{idx}.Cd;
            
            alpha_high = polar{idx+1}.Alpha;
            Cl_high_data = polar{idx+1}.Cl;
            Cd_high_data = polar{idx+1}.Cd;
            
            Cl_low = interp1(alpha_low, Cl_low_data, alpha_deg, 'linear', 'extrap');
            Cd_low = interp1(alpha_low, Cd_low_data, alpha_deg, 'linear', 'extrap');
            
            Cl_high = interp1(alpha_high, Cl_high_data, alpha_deg, 'linear', 'extrap');
            Cd_high = interp1(alpha_high, Cd_high_data, alpha_deg, 'linear', 'extrap');
            
            Cl = (1-w)*Cl_low + w*Cl_high;
            Cd = (1-w)*Cd_low + w*Cd_high;
            
            Cn = Cl*cos(phi) + Cd*sin(phi);
            Ct = Cl*sin(phi) - Cd*cos(phi);
            
            f = (B/2)*(R - r)/(r*sin(phi));
            F = (2/pi)*acos(exp(-f));
            F = max(F, 1e-4);
            
            a_new = 1 / ( (4*F*sin(phi)^2)/(sigma*Cn) + 1 );
            ap_new = 1 / ( (4*F*sin(phi)*cos(phi))/(sigma*Ct) - 1 );
            
            k = 0.2;
            a_j = (1-k)*a_j + k*a_new;
            ap_j = (1-k)*ap_j + k*ap_new;
            
            if abs(a_new - a_j) < tol && abs(ap_new - ap_j) < tol
                break;
            end
        end
        
        a(j) = a_j;
        ap(j) = ap_j;
        phi_arr(j) = rad2deg(phi);
        aoa_arr(j) = alpha_deg;
        
        dT(j) = 4*pi*rho*U^2*a_j*(1-a_j)*r*dr;
    end

    dQ = 4*pi*rho*U*(1-a).*ap*Omega.*r_vec.^3*dr;
    P = sum(Omega * dQ);

    A = pi*R^2;
    Cp = P / (0.5 * rho * A * U^3);

    Cp_list(L) = Cp;

    % Store best
    if Cp > best_Cp
        best_Cp = Cp;
        best_lambda = lambda;

        best_results.a = a;
        best_results.ap = ap;
        best_results.aoa_arr = aoa_arr;
        best_results.phi_arr = phi_arr;
        best_results.dT = dT;
    end
end

% Cp vs Lambda
figure;
plot(lambda_range, Cp_list, 'LineWidth',2);
hold on;
plot(lambda_range, Cp_list, '*');
xlabel('\lambda'); ylabel('C_p');
title('Cp vs Tip Speed Ratio');
grid on;

fprintf('\nOptimal lambda = %.2f\n', best_lambda);
fprintf('Maximum Cp = %.4f\n', best_Cp);

% Plots for best lambda

figure;
plot(r_vec, best_results.a, 'LineWidth',2); hold on;
plot(r_vec, best_results.a, '*');
xlabel('Radius (m)'); ylabel('a');
title('Axial Induction Factor vs Radius');
grid on;

figure;
plot(r_vec, best_results.ap, 'LineWidth',2); hold on;
plot(r_vec, best_results.ap, '*');
xlabel('Radius (m)'); ylabel('a''');
title('Tangential Induction Factor vs Radius'); grid on;

figure;
plot(r_vec, best_results.aoa_arr, 'LineWidth',2); hold on;
plot(r_vec, best_results.aoa_arr, '*');
xlabel('Radius (m)'); ylabel('AoA (deg)');
title('AOA Distribution vs Radius'); grid on;

figure;
plot(r_vec, best_results.phi_arr, 'LineWidth',2); hold on;
plot(r_vec, best_results.phi_arr, '*');
xlabel('Radius (m)'); ylabel('\phi (deg)');
title('Flow Angle vs Radius'); grid on;

figure;
plot(r_vec, best_results.dT, 'LineWidth',2); hold on;
plot(r_vec, best_results.dT, '*');
xlabel('Radius (m)'); ylabel('dT');
title('Thrust Distribution vs Radius'); grid on;