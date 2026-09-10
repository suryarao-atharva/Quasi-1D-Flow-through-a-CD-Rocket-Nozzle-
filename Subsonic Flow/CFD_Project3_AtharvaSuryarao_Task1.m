%Subsonic Flow

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
CFL = 0.25;

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
dlnA_dx = gradient(log(A), (dx));  % d/dx lnA

%Boundary Conditions at inlet and outlet
P_inlet = 6894.76;
P_back = 6136.34;
T_inlet = 55.56;

% Guessing the Initial Fields;
P = linspace(P_inlet,P_back,N);  % Initial guess of Pressure
P(1) = P_inlet;
P(N) = P_back;

V = linspace (20,30,N);  % Initial guess of velocity
T = T_inlet - (V.^2 / (2 * cp));  % Initial guess of density

rho = P./(R*T); % Initial guess of density

a = (gamma*R*T).^0.5; % Speed of Sound

iterations = 0; % iteration counter

% Initializing the Conservative variables provided in the project
Q1 = rho.*A;                                   
Q2 = rho.*A.*V;   
Q3 = rho.*A.*(cv.*T + 0.5.*V.^2);

% Initializing the fluxes provided in the project
F1 = Q2;
F2 = (Q2.^2./Q1) + (gamma-1) .* (Q3 - 0.5.*Q2.^2./ Q1);
F3 = (gamma.*Q2.*Q3./Q1) - 0.5.*(gamma - 1).*(Q2.^3./Q1.^2);

J1 = 0;
J2 = (gamma - 1).*(Q3 - 0.5.*Q2.^2./Q1).*dlnA_dx;
J3 = 0;

Q2_old_residual = Q2;

while residual > tolerance

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

    % Recover temperature, density, velocity, Pressure for interior nodes
    % and then apply the boundary conditions for nodes 1 and N

    rho = Q1./A;
    V = Q2./Q1;
    T = (Q3./Q1 - 0.5.* V.^2)./cv;
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

    P(N) = P_back;
    V(N) = 2*V(N-1)- V(N-2); % Extrapolating the velocity at outlet
    T(N) = 2*T(N-1)- T(N-2); % Extrapolating the temperature at outlet
    rho(N) = P(N)/(R*T(N));

    % Update values of Q1,Q2,Q3 at the outlet using the updated values of primitive variables at outlet
    Q1(N) = rho(N) * A(N);                                   
    Q2(N) = rho(N) * A(N) * V(N);   
    Q3(N) = rho(N) * A(N) * (cv*T(N) + 0.5*V(N)^2);

    % Calculating the residual
    residual = norm(Q2 - Q2_old_residual) / (norm(Q2_old_residual) + eps); % Normalized Norm
    Q2_old_residual = Q2;
end

% Getting the final results

Pressure_ratio = P ./P_inlet;
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
    0.2305986, 0.2311824, 0.2329429, 0.2359077, 0.2401235, ...
    0.2456567, 0.2525952, 0.2610493, 0.2711545, 0.2830723, ...
    0.2969925, 0.3131352, 0.3317496, 0.3531143, 0.3775311, ...
    0.4053157, 0.4367796, 0.4721986, 0.5117609, 0.5554809, ...
    0.6030496, 0.6535671, 0.7050224, 0.7532675, 0.7904066, ...
    0.8049845, 0.7999743, 0.7858996, 0.7650277, 0.7398368, ...
    0.7123653, 0.6840839, 0.6559903, 0.6287439, 0.6027672, ...
    0.5783231, 0.5555651, 0.5345717, 0.5153693, 0.4979494, ...
    0.4822819, 0.4683219, 0.4560181, 0.4453161, 0.4361631, ...
    0.4285094, 0.4223110, 0.4175303, 0.4141374, 0.4121103, ...
    0.4114360];

% Pressure ratio

p_Pa = [ ...
    6.644140e+03, 6.642900e+03, 6.639142e+03, 6.632755e+03, 6.623548e+03, ...
    6.611241e+03, 6.595458e+03, 6.575702e+03, 6.551338e+03, 6.521571e+03, ...
    6.485420e+03, 6.441674e+03, 6.388867e+03, 6.325248e+03, 6.248782e+03, ...
    6.157150e+03, 6.047850e+03, 5.918374e+03, 5.766540e+03, 5.591082e+03, ...
    5.392589e+03, 5.175028e+03, 4.948327e+03, 4.732902e+03, 4.566198e+03, ...
    4.500734e+03, 4.523227e+03, 4.586443e+03, 4.680153e+03, 4.793054e+03, ...
    4.915678e+03, 5.041067e+03, 5.164454e+03, 5.282679e+03, 5.393787e+03, ...
    5.496651e+03, 5.590737e+03, 5.675919e+03, 5.752347e+03, 5.820352e+03, ...
    5.880361e+03, 5.932857e+03, 5.978327e+03, 6.017243e+03, 6.050045e+03, ...
    6.077122e+03, 6.098812e+03, 6.115391e+03, 6.127078e+03, 6.134028e+03, ...
    6.136336e+03 ...
];


% Mach Number for analytical 
M_anyl = [ ...
0.2305986, 0.2311824, 0.2329429, 0.2359077, 0.2401235, ...
0.2456567, 0.2525952, 0.2610493, 0.2711545, 0.2830723, ...
0.2969925, 0.3131352, 0.3317496, 0.3531143, 0.3775311, ...
0.4053157, 0.4367796, 0.4721986, 0.5117609, 0.5554809, ...
0.6030496, 0.6535671, 0.7050224, 0.7532675, 0.7904066, ...
0.8049845, 0.7999743, 0.7858996, 0.7650277, 0.7398368, ...
0.7123653, 0.6840839, 0.6559903, 0.6287439, 0.6027672, ...
0.5783231, 0.5555651, 0.5345717, 0.5153693, 0.4979494, ...
0.4822819, 0.4683219, 0.4560181, 0.4453161, 0.4361631, ...
0.4285094, 0.4223110, 0.4175303, 0.4141374, 0.4121103, ...
0.4114360];

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
plot(x_m, P_anyl, 'b-o', 'LineWidth', 1.5);
plot(x_m, Pressure_ratio, 'r', 'LineWidth', 1.5);
xlabel('x (m)');
ylabel("Pressure")
legend("Analytical Solution", "Code Data")
title('Comparison of Pressure with data with Analytical Solution');
grid on;

%Density ratio Comparison
figure()
hold on;
plot(x_m, rho_anyl, 'b-o', 'LineWidth', 1.5);
plot(x_m, rho_ratio, 'r', 'LineWidth', 1.5);
xlabel('x (m)');
ylabel("density ratio")
legend("Analytical Solution", "Code Data")
title('Comparison of density with data with Analytical Solution');
grid on;

%Temperature ratio Comparison
figure()
hold on;
plot(x_m, T_anyl, 'b-o', 'LineWidth', 1.5);
plot(x_m, Temperature_ratio, 'r', 'LineWidth', 1.5);
xlabel('x (m)');
ylabel("Temperature")
legend("Analytical Solution", "Code Data")
title('Comparison of Temperature with data with Analytical Solution');
grid on;

%Mach Number ratio Comparison
figure()
hold on;
plot(x_m, M_anyl, 'b-o', 'LineWidth', 1.5);
plot(x_m, Mach, 'r', 'LineWidth', 1.5);
xlabel('x (m)');
ylabel("Mach number")
legend("Analytical Solution", "Code Data")
title('Comparison of Mach number with data with Analytical Solution');
grid on;


% Comparison of Pressure with NASA's Data
figure()
hold on;
plot(x_m, p_Pa./P_inlet, 'b-o', 'LineWidth', 1.5);
plot(x_m, P./P_inlet, 'r', 'LineWidth', 1.5);
xlabel('x (m)');
ylabel("Pressure")
legend("NASA", "Code Data")
title('Comparison of Pressure with data from NASA');
grid on;

% Comparison of Mach Number with NASA's Data
figure()
hold on;
plot(x_m, M_NASA, 'b-o', 'LineWidth', 1.5);
plot(x_m, Mach, 'r', 'LineWidth', 1.5);
xlabel('x (m)');
ylabel("Mach Number")
legend("NASA", "Code Data")
title('Comparison of Mach Number with data from NASA');
grid on;