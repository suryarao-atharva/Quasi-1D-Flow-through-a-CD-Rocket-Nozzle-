clc;
clear all;

% Simulation of fluid flow through a CD nozzle for a supersonic flow using MacCormack's Scheme 

% Define parameters 
n = 51; % number of points in the domain
L = 0.254; % Length of the Nozzle
delta_x = L/(n-1); % Grid Spacing
x = linspace(0, L, n); % gives the x-coordinate for the points in the nozzle

% disp(x) % we can see the x coordinates of our nozzle

% Gas properties
gamma = 1.4; % heat capacity ratio
R = 287; % units are in J/kg.k
Cv = R/(gamma-1); 
CFL = 0.2; % Courant Number required for stability
tolerance = 1e-5; % to check the convergence
max_iterations = 1e+5;


% Values of Pressure and Temperature are known at the inlet and Outlet

P_inlet = 6894.76; % Units are in Pa
P_outlet = 1103.16; % Units are in Pa
T_inlet = 55.56; % Units are in K
rho_inlet = P_inlet/(R*T_inlet); % Density at inlet is constant

% Adding the Area function
A = zeros(1, n);   % Area
% 
for i = 1:n
    if x(i)>=0 && x(i)< 0.127
        A(i) = 0.0444 - 0.019*cos((pi)*((0.2*x(i)/0.0254)- 1));

    elseif x(i)>= 0.127 && x(i)<=0.254
        A(i) = 0.0318 - 0.0063*cos((pi)*((0.2*x(i)/0.0254)- 1));
    end
end
% disp(A)


% Initialization of dimensionalized nozzle flow variables
rho = zeros(1,n);
T = zeros(1,n);
V = zeros(1,n);
P = zeros(1,n);


Q2_inlet_guess = 0.7;

% Setting up the initial conditions
for i  = 1:n
    rho(i) = rho_inlet*(1 - (0.5/L)*x(i)); %Initial guess of Density in (kg/m^3)
    T(i) = T_inlet*(1 - (0.5/L)*x(i)); %Initial guess of Temperature in (K)
    V(i) = (Q2_inlet_guess)/(rho(i) * A(i)); %Initial guess of velocity
    P(i) = rho(i)*R*T(i);
end

% The inlet pressure and outlet pressure for the subsonic case is mentioned and hence in the 
% initial Pressure guess we can set the Pressure at i = 1, n as the the inlet and back Pressure 
P(1) = P_inlet;
% P(n) = P_outlet;

% Initialization of Conservative Variables

Q1 = rho.* A; % we get the initial array for Q1
Q2 = rho.* A.* V; % we get the initial array for Q2
Q3 = rho.* A.*(Cv*T + 0.5*V.^2); % we get the initial array for Q2

%Initialization of Flux Variables

F1 = Q2;
F2 = ((Q2.^2)./(Q1)) + (gamma-1) * (Q3 - 0.5*((Q2.^2)./(Q1)));
F3 = (Q2./Q1).*(Q3 + (gamma-1) * (Q3 - 0.5*((Q2.^2)./(Q1))));

%Initialization of Source Terms
J1 = 0;
J3 = 0;

% derivative d(ln A)/dx (needed for source J2)
dlnA_dx = gradient(log(A), (delta_x));  % d/dx ln A

J2 = ((gamma-1) * (Q3 - 0.5*((Q2.^2)./(Q1)))).*(dlnA_dx);

% Initialization for the McCormacks's scheme Predictor and Corrector terms

% Initialization of Predictor step derivative terms
dQ1_dt_predicted = zeros(1, n);
dQ2_dt_predicted = zeros(1, n);
dQ3_dt_predicted = zeros(1, n);

% Initialization of Predictor step derivative terms
Q1_predicted = zeros(1, n);
Q2_predicted = zeros(1, n);
Q3_predicted = zeros(1, n); 

% Initialization of Corrector step derivative terms
dQ1_dt_corrected = zeros(1, n);
dQ2_dt_corrected = zeros(1, n);
dQ3_dt_corrected = zeros(1, n);


residual = 10;
iterations = 0;

while residual > tolerance
    iterations = iterations+1;

    %Computing time step for CFL Condition
    a = sqrt(gamma * R * T);  % local speed of sound
    delta_t = 0.9*(CFL * delta_x) / (max(abs(V) + a));

    
    % Predictor Step (FTFS)
    for j = 1:n-1
    
        dQ1_dt_predicted(j) = -(F1(j+1)-F1(j))/delta_x;
        dQ2_dt_predicted(j) = J2(j) -(F2(j+1)-F2(j))/delta_x;
        dQ3_dt_predicted(j) = -(F3(j+1)-F3(j))/delta_x;
        
        Q1_predicted(j) = Q1(j) + delta_t * dQ1_dt_predicted(j);
        Q2_predicted(j) = Q2(j) + delta_t * dQ2_dt_predicted(j);
        Q3_predicted(j) = Q3(j) + delta_t * dQ3_dt_predicted(j);
        
    end

    % Q1_predicted(1) = Q1(1);
    % Q2_predicted(1) = Q2(1);
    % Q3_predicted(1) = Q3(1);
    
    % Q1_predicted(1) = rho_inlet * A(1);
    % Q2_predicted(1) = 2*Q2_predicted(2) - Q2_predicted(3);
    % 
    % V_predict(1) = Q2_predicted(1)/Q1_predicted(1);
    % Q3_predicted(1) = rho_inlet * A(1) * (Cv*T_inlet + 0.5*V_predict(1)^2);


    Q1_predicted(n) = 2*Q1_predicted(n-1)- Q1_predicted(n-2);
    Q2_predicted(n) = 2*Q2_predicted(n-1)- Q2_predicted(n-2);
    Q3_predicted(n) = 2*Q3_predicted(n-1)- Q3_predicted(n-2);
    
    % For the Corrector step we need the updated value of Flux and Source terms which are also predicted values
    
    F1_predicted = Q2_predicted;
    F2_predicted = ((Q2_predicted.^2)./(Q1_predicted)) + (gamma-1) * (Q3_predicted - 0.5*((Q2_predicted.^2)./(Q1_predicted)));
    F3_predicted = (Q2_predicted./Q1_predicted).*(Q3_predicted + (gamma-1) * (Q3_predicted - 0.5*((Q2_predicted.^2)./(Q1_predicted)))); 
    
    J2_predicted = ((gamma-1) * (Q3_predicted - 0.5*((Q2_predicted.^2)./(Q1_predicted)))).*(dlnA_dx);
    
    % Corrector step (FTBS)
    
    for k = 2:n-1
        dQ1_dt_corrected(k) = - (F1_predicted(k)-F1_predicted(k-1))/delta_x;
        dQ2_dt_corrected(k) = J2_predicted(k) - (F2_predicted(k)-F2_predicted(k-1))/delta_x;
        dQ3_dt_corrected(k) = - (F3_predicted(k)-F3_predicted(k-1))/delta_x;
    end
    
    % Average of the derivative of conservative variables
    
    average_dQ1_dt = 0.5.*(dQ1_dt_predicted + dQ1_dt_corrected);
    average_dQ2_dt = 0.5.*(dQ2_dt_predicted + dQ2_dt_corrected);
    average_dQ3_dt = 0.5.*(dQ3_dt_predicted + dQ3_dt_corrected);
    
    % updated value after a time step of the conservative variables
    
    Q1_new = Q1 + delta_t * average_dQ1_dt;
    Q2_new = Q2 + delta_t * average_dQ2_dt;
    Q3_new = Q3 + delta_t * average_dQ3_dt;
    
    % disp("Q3_new")
    % disp(Q3_new)
    % Now we have the Correct values of Q1,Q2,Q3 but we still need to correct the values at the boundaries for a subsonic flow problem
    
    % Inlet Boundary Conditions
    Q1_new(1) = rho_inlet*A(1);
    Q2_new(1) = 2*Q2_new(2) - Q2_new(3);
    
    V_new(1) = Q2_new(1)/Q1_new(1);
    Q3_new(1) = rho_inlet * A(1) * (Cv*T_inlet + 0.5*V_new(1)^2);

    % Q3_new(1) = (P_inlet*A(1))/(gamma-1) + 0.5*Q2_new(1)^2/Q1_new(1);
    
    % Outlet Boundary Condtions
    
    Q1_new(n) = 2*Q1_new(n-1) - Q1_new(n-2);
    Q2_new(n) = 2*Q2_new(n-1) - Q2_new(n-2);
    Q3_new(n) = 2*Q3_new(n-1) - Q3_new(n-2);

    
    % disp("Q3_new after updating boundary")
    % disp(Q3_new)
    
    %Updating the Flux terms
    
    % First we update the Flux terms
    F1 = Q2_new;
    F2 = ((Q2_new.^2)./(Q1_new)) + (gamma-1) * (Q3_new - 0.5*((Q2_new.^2)./(Q1_new)));
    F3 = (Q2_new./Q1_new).*(Q3_new + (gamma-1) * (Q3_new - 0.5*((Q2_new.^2)./(Q1_new))));
    
    % Next we update the Source terms
    J2 = ((gamma-1) * (Q3_new - 0.5*((Q2_new.^2)./(Q1_new)))).*(dlnA_dx);
    
    % new primitives variables calculation
    rho = Q1_new ./ A;
    V = Q2_new ./ Q1_new;
    T = ((Q3_new ./ (rho .* A)) - 0.5 .* V.^2) ./ Cv;
    P = rho .* R .* T;
    M = V ./ sqrt(gamma * R .* T);
    
    % convergence check 

    residual_1 = max(abs(Q1_new - Q1) ./ (abs(Q1_new) + eps));
    residual_2 = max(abs(Q2_new - Q2) ./ (abs(Q2_new) + eps));
    residual_3 = max(abs(Q3_new - Q3) ./ (abs(Q3_new) + eps));
    residual = max([residual_1, residual_2, residual_3]);

    % Updating Q
    Q1 = Q1_new;
    Q2 = Q2_new;
    Q3 = Q3_new;
end

% Plot results 

figure();
plot(x, P, 'LineWidth', 1.5);
title('Pressure Distribution along the Nozzle');
xlabel('Nozzle X-Distance');
ylabel(' Pressure');
grid on;

figure();
plot(x, rho, 'LineWidth', 1.5);
title('Density Distribution along the Nozzle');
xlabel('Nozzle X-Distance');
ylabel(' Density');
grid on;

figure();
plot(x, T, 'LineWidth', 1.5);
title('Temperature Distribution along the Nozzle');
xlabel('Nozzle X-Distance');
ylabel('Non-dimensional Temperature');
grid on;





x_nozzle_nasa = [0.0000 0.00508 0.01016 0.01524 0.02032 0.02540 0.03048 0.03556 0.04064 ...
     0.04572 0.05080 0.05588 0.06096 0.06604 0.07112 0.07620 0.08128 0.08636 ...
     0.09144 0.09652 0.10160 0.10668 0.11176 0.11684 0.12192 0.12700 0.13208 ...
     0.13716 0.14224 0.14732 0.15240 0.15748 0.16256 0.16764 0.17272 0.17780 ...
     0.18288 0.18796 0.19304 0.19812 0.20320 0.20828 0.21336 0.21844 0.22352 ...
     0.22860 0.23368 0.23876 0.24384 0.24892 0.25400];

M_nasa = [0.2395428 0.2401525 0.2419912 0.2450882 0.2494931 0.2552769 0.2625332 ...
     0.2713805 0.2819642 0.2944594 0.3090735 0.3260485 0.3456640 0.3682380 ...
     0.3941266 0.4237210 0.4574398 0.4957166 0.5389794 0.5876219 0.6419665 ...
     0.7022194 0.7684249 0.8404243 0.9178230 0.9997698 1.049278 1.099540 ...
     1.150450 1.201643 1.252754 1.303419 1.353288 1.402030 1.449334 1.494914 ...
     1.538515 1.579906 1.618885 1.655275 1.688924 1.719701 1.747494 1.772210 ...
     1.793771 1.812114 1.827187 1.838950 1.847373 1.852435 1.854124];

figure();
hold on
plot(x_nozzle_nasa, M_nasa,"LineStyle",":","LineWidth",2,"Color",'red');
plot(x, M, 'LineWidth', 1.5,"Color",'blue',"LineWidth",2);

title('Mach Number along the Nozzle');
xlabel('Nozzle X-Distance');
ylabel('Mach Number');
grid on;



