% Supersonic Flow
clc;
clear all;

N = 51; % Number of grid points
L = 0.254; % Length of the Nozzle
x = linspace(0,L,N); % Gives us the x coordinates
dx = L/(N-1); % Spacing between two grid points

%Gas Properties
R = 287;
gamma = 1.4;
cv = R/(gamma-1);
cp = gamma*cv;

tolerance = 1e-6;
residual = 10;
CFL = 0.5;

%Area function
A = zeros(1,N);
for i = 1:N
    if x(i)<0.127
        A(i) = 0.0444 - 0.019 * cos((0.2 * x(i) / 0.0254 - 1) * pi);
    else
        A(i) = 0.0318 - 0.0063 * cos((0.2 * x(i) / 0.0254 - 1) * pi);

    end
end

%Derivative of the Area function (Required for J2)
dlnA_dx = gradient(log(A), (dx));  % d/dx ln A

%Boundary Conditions at inlet and outlet
P_inlet = 6894.76;
P_back = 1103.16;
T_inlet = 55.56;

% Guessing the Initial Fields;
P = linspace(P_inlet,P_back,N);  % Initial guess of Pressure
V = linspace (20,300,N);  % Initial guess of velocity increases rapidly for supersonic flow

T = T_inlet - (V.^2 / (2 * cp));  % Initial guess of density
T(T <= 0) = T_inlet*0.5; % To ensure no temperature is lesser than zero

rho = P./(R*T); % Initial guess of density
a = (gamma*R*T).^0.5; % Speed of Sound

iterations = 0;
max_iterations = 1e+5;

% Initializing the Conservative variables provided in the project
Q1 = rho.*A;                                   
Q2 = rho.*A.*V;   
Q3 = rho.*A.*(cv.*T + 0.5.*V.^2);

% Initializing the fluxes provided in the project
F1 = Q2;
F2 = (Q2.^2./Q1) + (gamma-1) .* (Q3 - 0.5.*Q2.^2./ Q1);
F3 = (gamma.*Q2.*Q3./Q1) - 0.5.*(gamma - 1).*(Q2.^3./Q1.^2);

J1 = zeros(1,N);
J2 = (gamma - 1).*(Q3 - 0.5.*Q2.^2./Q1).*dlnA_dx;
J3 = zeros(1,N);

Q2_old_residual = Q2;

while residual > tolerance
    iterations  = iterations+1;

    Q1_old = Q1;
    Q2_old = Q2;
    Q3_old = Q3;

    % Dynamic time step (expression derived from John Anderson Book)
    dt = (CFL*dx)/(max(abs(a)+ abs(V)));

    F1 = Q2;
    F2 = (Q2.^2./Q1) + (gamma-1).*(Q3 - 0.5.*Q2.^2./ Q1);
    F3 = (gamma.*Q2.*Q3./Q1) - 0.5.*(gamma - 1).*(Q2.^3./Q1.^2);
    
    J2 = (gamma - 1).*(Q3 - 0.5.*Q2.^2./Q1).*dlnA_dx;

    % Predictor Step (1st step in MacCormack's Scheme) (Forward Time Forward Space discretization)

    Q1_star = Q1;
    Q2_star = Q2;
    Q3_star = Q3;

    for j = 2:N-1
        Q1_star(j) = Q1(j) - (dt/dx) * (F1(j+1) - F1(j)); % J1 term is 0
        Q2_star(j) = Q2(j) - (dt/dx) * (F2(j+1) - F2(j)) + dt * J2(j);
        Q3_star(j) = Q3(j) - (dt/dx) * (F3(j+1) - F3(j)); % J3 term is 0
    end

    % Now we update the values of Fluxes and Source term using the
    % predicted values of Q1_star, Q2_star, Q3_star

    F1_star = Q2_star;
    F2_star = (Q2_star.^2./Q1_star) + (gamma-1).*(Q3_star - 0.5.*Q2_star.^2./ Q1_star);
    F3_star = (gamma.*Q2_star.*Q3_star./Q1_star) - 0.5.*(gamma - 1).*(Q2_star.^3./Q1_star.^2);
    
    J2_star = (gamma - 1).*(Q3_star - 0.5.*Q2_star.^2./Q1_star).*dlnA_dx;

    % Corrector Step (2nd step in MacCormack's Scheme) (Forward Time Backward Space discretization)
    for k = 2:N-1
        Q1(k) = 0.5 * (Q1_old(k) + Q1_star(k) - (dt/dx)*(F1_star(k) - F1_star(k-1)));
        Q2(k) = 0.5 * (Q2_old(k) + Q2_star(k) - (dt/dx)*(F2_star(k) - F2_star(k-1)) + dt * J2_star(k));
        Q3(k) = 0.5 * (Q3_old(k) + Q3_star(k) - (dt/dx)*(F3_star(k) - F3_star(k-1)));
    end

    % Recover temperature, density, velocity, Pressure for interior nodes and then apply the boundary conditions for nodes 1 and N
    rho = Q1./A;
    V = Q2./Q1;
    T = (Q3./Q1 - 0.5.*V.^2)./cv;
    P = rho.*R.*T;
    a = sqrt(gamma*R.*T);

    % Applying the boundary Conditions at inlet

    V(1) = 2*V(2) - V(3);  % Extrapolating the velocity at inlet
    T(1) = T_inlet - V(1)^2 / (2*cp);
    P(1) = P_inlet * (T(1)/T_inlet)^(gamma/(gamma - 1));
    rho(1) = P(1) / (R*T(1));

    % Update values of Q1,Q2,Q3 at the inlet using the updated values of primitive variables at inlet
    Q1(1) = rho(1) * A(1);                                   
    Q2(1) = rho(1) * A(1) * V(1);   
    Q3(1) = rho(1) * A(1) * (cv*T(1) + 0.5*V(1)^2);
    
    % Applying the boundary Conditions at outlet
    rho(N) = 2*rho(N-1)- rho(N-2); % Extrapolating the density at outlet
    V(N) = 2*V(N-1)- V(N-2); % Extrapolating the velocity at outlet
    T(N) = 2*T(N-1)- T(N-2); % Extrapolating the temperature at outlet

    
    %getting the pressure at node N

    P(N) = rho(N)*R*T(N);

    % Update values of Q1,Q2,Q3 at the outlet using the updated values of primitive variables at outlet
    Q1(N) = rho(N) * A(N);                                   
    Q2(N) = rho(N) * A(N) * V(N);   
    Q3(N) = rho(N) * A(N) * (cv*T(N) + 0.5*V(N)^2);

    % Calculating the residual

    residual = norm(Q2 - Q2_old_residual) / (norm(Q2_old_residual) + eps);
    Q2_old_residual = Q2;
end

% Getting the final results
Pressure_ratio = P ./ P_inlet;
rho_ratio = rho ./ (P_inlet / (R * T_inlet));
Temperature_ratio = T ./ T_inlet;
Mach = V ./ a;

% NASA's data for Pressure ratio for a subsonic Flow;

% axial Coordinates for the nozzle
x_m = [0.000000, 0.005080, 0.010160, 0.015240, 0.020320, 0.025400, 0.030480, 0.035560, 0.040640, 0.045720,...
    0.050800, 0.055880, 0.060960, 0.066040, 0.071120, 0.076200, 0.081280, 0.086360, 0.091440, 0.096520,...
    0.101600, 0.106680, 0.111760, 0.116840, 0.121920, 0.127000, 0.132080, 0.137160, 0.142240, 0.147320,...
    0.152400, 0.157480, 0.162560, 0.167640, 0.172720, 0.177800, 0.182880, 0.187960, 0.193040, 0.198120,...
    0.203200, 0.208280, 0.213360, 0.218440, 0.223520, 0.228600, 0.233680, 0.238760, 0.243840, 0.248920, 0.254000];


% Mach Number for NASA
M_NASA = [ ...
    0.2395428, 0.2401525, 0.2419912, 0.2450882, 0.2494931, ...
    0.2552769, 0.2625332, 0.2713805, 0.2819642, 0.2944594, ...
    0.3090735, 0.3260485, 0.3456640, 0.3682380, 0.3941266, ...
    0.4237210, 0.4574398, 0.4957166, 0.5389794, 0.5876219, ...
    0.6419665, 0.7022194, 0.7684249, 0.8404243, 0.9178230, ...
    0.9997698, 1.049278, 1.099540, 1.150450, 1.201643, ...
    1.252754, 1.303419, 1.353288, 1.402030, 1.449334, ...
    1.494914, 1.538515, 1.579906, 1.618885, 1.655275, ...
    1.688924, 1.719701, 1.747494, 1.772210, 1.793771, ...
    1.812114, 1.827187, 1.838950, 1.847373, 1.852435, ...
    1.854124 ...
    ];


% Pressure ratio
p_Pa = [ ...
    6.624825e+03, 6.623483e+03, 6.619422e+03, 6.612517e+03, 6.602563e+03, ...
    6.589253e+03, 6.572173e+03, 6.550782e+03, 6.524387e+03, 6.492109e+03, ...
    6.452862e+03, 6.405303e+03, 6.347792e+03, 6.278345e+03, 6.194623e+03, ...
    6.093898e+03, 5.963318e+03, 5.810302e+03, 5.632164e+03, 5.426404e+03, ...
    5.190767e+03, 4.923308e+03, 4.622400e+03, 4.296748e+03, 3.945457e+03, ...
    3.644037e+03, 3.437446e+03, 3.229693e+03, 3.030303e+03, 2.837864e+03, ...
    2.651063e+03, 2.472645e+03, 2.302326e+03, 2.140838e+03, 1.988911e+03, ...
    1.847261e+03, 1.716589e+03, 1.597569e+03, 1.489844e+03, 1.392916e+03, ...
    1.306073e+03, 1.228480e+03, 1.159243e+03, 1.097363e+03, 1.041771e+03, ...
    9.913823e+02, 9.451218e+02, 9.011752e+02, 8.619175e+02, 8.255921e+02, ...
    7.916280e+02 ...
];


% Mach Number for analytical
M_anyl = [ ...
0.2395428, 0.2401525, 0.2419912, 0.2450882, 0.2494931, ...
0.2552769, 0.2625332, 0.2713805, 0.2819642, 0.2944594, ...
0.3090735, 0.3260485, 0.3456640, 0.3682380, 0.3941266, ...
0.4237210, 0.4574398, 0.4957166, 0.5389794, 0.5876219, ...
0.6419665, 0.7022194, 0.7684249, 0.8404243, 0.9178230, ...
0.9997698, 1.049278, 1.099540, 1.150450, 1.201643, ...
1.252754, 1.303419, 1.353288, 1.402030, 1.449334, ...
1.494914, 1.538515, 1.579906, 1.618885, 1.655275, ...
1.688924, 1.719701, 1.747494, 1.772210, 1.793771, ...
1.812114, 1.827187, 1.838950, 1.847373, 1.852435, ...
1.854124 ...
];

% analytical solution for pressure
P_anyl = (1 + 0.2* M_anyl.^2).^-(1.4/0.4);

% analytical solution for temperature
T_anyl = 1./(1 + 0.2*M_anyl.^2);

% analytical solution for density
rho_anyl = (1 + 0.2* M_anyl.^2).^-(1/0.4);

% Comparison of results with Analytical Solution by plotting them

%Pressure ratio Comparison
figure()
hold on;
plot(x_m, P_anyl, 'g-o', 'LineWidth', 1.5);
plot(x_m, Pressure_ratio, 'r', 'LineWidth', 2);
xlabel('x (m)');
ylabel("Pressure")
legend("Analytical Solution", "Code Data")
title('Comparison of Pressure with data with Analytical Solution');
grid on;

%Density ratio Comparison
figure()
hold on;
plot(x_m, rho_anyl, 'g-o', 'LineWidth', 1.5);
plot(x_m, rho_ratio, 'r', 'LineWidth', 2);
xlabel('x (m)');
ylabel("density ratio")
legend("Analytical Solution", "Code Data")
title('Comparison of density with data with Analytical Solution');
grid on;

%Temperaure ratio Comparison
figure()
hold on;
plot(x_m, T_anyl, 'g-o', 'LineWidth', 1.5);
plot(x_m, Temperature_ratio, 'r', 'LineWidth', 2);
xlabel('x (m)');
ylabel("Temperature")
legend("Analytical Solution", "Code Data")
title('Comparison of Temperature with data with Analytical Solution');
grid on;

%Mach Number ratio Comparison
figure()
hold on;
plot(x_m, M_anyl, 'g-o', 'LineWidth', 1.5);
plot(x_m, Mach, 'r', 'LineWidth', 2);
xlabel('x (m)');
ylabel("Mach number")
legend("Analytical Solution", "Code Data")
title('Comparison of Mach number with data with Analytical Solution');
grid on;


% Comparison of Pressure with NASA's Data
figure()
hold on;
plot(x_m, p_Pa./P_inlet, 'g-o', 'LineWidth', 1.5);
plot(x_m, P./P_inlet, 'r', 'LineWidth', 2);
xlabel('x (m)');
ylabel("Pressure")
legend("NASA", "Code Data")
title('Comparison of Pressure with data from NASA');
grid on;

% Comparison of Mach Number with NASA's Data
figure()
hold on;
plot(x_m, M_NASA, 'g-o', 'LineWidth', 1.5);
plot(x_m, Mach, 'r', 'LineWidth', 2);
xlabel('x (m)');
ylabel("Mach Number")
legend("NASA", "Code Data")
title('Comparison of Mach Number with data from NASA');
grid on;