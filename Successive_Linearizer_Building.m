function [A, B, C, D, U, Y, X, DX, poles] = Successive_Linearizer_Building(x, u, d)
%#codegen

% Define constant outputs
Ts = 0.5;
C = eye(2);
D = zeros(2,2);
% Nominal U are obtained from measurements
U = [u1; u2];
% Nominal X and Y are obtained from estimated MPC states
Y = x;
X = y;
% Analytical linearization of mechanistic CSTR model (continuous time)
[A, Bo] = getContinuous(x, u, d);
% Convert continuous to discrete time
[A, Bo] = getDiscrete(A, Bo, Ts);
% Last column of B contains the constant offset contribution.
DX = Bo(:,3);
B = Bo(:,1:2);
% Compute poles
poles = abs(eig(A));

function [a, b] = getContinuous(X, U, Dist)
% Define linear continuous-time matrices A and B for CSTR.  Last input is
% offset.
K = zeros(15,1);
% Define constants
% Overall Model Constantsg = 9.81 * 1e2 * 3600^2;           % Gravitational constant [cm h^-2]
g = 9.81 * 1e2 * 3600^2;           % Gravitational constant [cm h^-2]
K_valve = 0.865;    % Valve coefficient conversion constant [to SI units]

rho_water = 1000 * 1e-6;   % Density of water [kg cm^-3]
rho_oil = 800 * 1e-6;      % Density of heating oil [kg cm^-3]
Cp_water = 4180 / 3600;    % Specific heat capacity of water [Wh kg^-1 K^-1]
Cp_oil = 1700 / 3600;      % Specific heat capacity of heating oil [Wh kg^-1 K^-1]

% Drain Water Heat Recovery (DWHR) Constants
n = 30;             % Number of internal tubes
d_tube = 0.01 * 1e2;      % Internal tube diameter [cm]
D_tube = 0.15 * 1e2;      % External tube diameter [cm]
L = 5 * 1e2;              % Tube(s) length [cm]
A_s_dwhr = n * pi * d_tube * L;              % Total surface area of internal tubes
V_dwhr_int = 0.25 * n * pi * L * d_tube^2;   % Total volume of internal tubes [m^3]
V_dwhr_tot = 0.25 * pi * L * D_tube^2;       % Volume of external tube [m^3]
V_dwhr_ext = V_dwhr_tot - V_dwhr_int;   % Volume of external tubes sans internal tubes [m^3]
V_dwhr_cold = V_dwhr_int;              % Cold water within internal tubes (for maintainability)
V_dwhr_hot = V_dwhr_ext;               % Hot water in external tube
U_dwhr = 1275 * 1e-4;                          % Overall heat transfer coefficient [W cm^-2 K^-1]

% Preheat Tank (PHT) Constants
H_pht_0 = 0.3 * 1e2;      % Initial level of PHT [cm]
H_pht_max = 1 * 1e2;      % Max height of PHT [cm]
A_pht = 0.01 * 1e4;       % Area of PHT [cm^2]
A_s_pht = 0.001 * 1e4;    % TODO Surface area of heating tube in PHT [cm^2]
U_pht = 1275 * 1e-4;       % TODO Overall heat transfer coefficient (PHT-LBTC) [W cm^-2 K^-1]
C_v_pht = 0.0167;   % Outlet valve coefficient

K(1) = U_dwhr * A_s_dwhr / (2 * rho_water * V_dwhr_cold * Cp_water); % [s^-1]
K(2) = K_valve * C_v_dhw * lambda / V_dwhr_hot; % [m^(-1/2) s^-1]
K(3) = Q_oil_tc / V_tube_tc; %[s^-1]
K(4) = 1 / (rho_oil * V_tube_tc * Cp_oil); % [J^-1 K]
K(5) = U_pht * A_s_pht / (rho_oil * V_tube_tc * Cp_oil); % [s^-1]
K(6) = K_valve * C_v_pht * lambda / A_pht; % [m^(-1/2) s^-1]
K(7) = U_pht * A_s_pht / (rho_water * Cp_water * A_pht); % [m s^-1]
K(8) = K_valve * C_v_pht * lambda / A_dhw; % [m^(-1/2) s^-1]
K(9) = K_valve * C_v_dhw * lambda / A_dhw; % [m^(-1/2) s^-1]
K(10) = 1 / (A_dhw * rho_water * Cp_water); % [J^-1 m K]
K(11) = Q_oil_tst / V_tube_tst; %[s^-1]
K(12) = U_tst * A_s_tst / (rho_water * V_tube_tst * Cp_oil); % [s^-1]
K(13) = 1 / (rho_water * Cp_water * V_tst); % [J^-1 K]
K(14) = 1 / V_dwhr_cold; % [m^-3]
K(15) = 1 / A_pht; % [m^-2]

% Nonlinear model time derivatives
dxdt = [(K(14)*U(1)+K(1))*Dist(1) - (K(14)*U(1)+K(1))*X(1)+ K(1)*X(8)-K(1)*X(2);
     -K(1)*Dist(1)+K(1)*X(1) + (K(2)*Dist(2)*sqrt(X(7)) - K(1))*X(8) - (K(2)*Dist(2)*sqrt(X(7)) - K(1))*X(2);
     K(15)*U(1) - K(6)*U(3)*sqrt(X(3));
     K(15)*U(1)*X(1)/X(3) - (K(15)*U(1)/X(3) + K(7)/X(3))*X(4) + K(7)*X(6)/X(3);
     K(3)*U(2)*X(6) - K(3)*U(2)*X(5) + K(4)*U(5);
     -(K(3)*U(2) + K(5))*X(6) + K(3)*U(2)*X(5) + K(5)*U(4);
     K(8)*U(3)*sqrt(X(3)) - K(9)*Dist(2)*sqrt(X(7));
     K(8)*U(3)*sqrt(X(3))*X(4)/X(7)-K(8)*U(3)*sqrt(X(3))*X(8)/X(7) + K(10)*U(7)/X(7);
     K(13)*U(6) - K(13)*U(7) + K(13) * U(8);
     K(11)*U(4)*X(11) - K(11)*U(4)*X(10) + K(4)*U(6);
     K(11)*U(4)*X(10) - (K(11)*U(4) + K(12))*X(11) + K(12)*X(9)];

a = [-K(1)-K(14)*U(1) -K(1) 0 0 0 0 0 K(1) 0 0 0;
    K(1) K(1)-Dist(2)*K(2)*sqrt(X(7)) 0 0 0 0 (Dist(2)*K(2)*X(8))/(2*sqrt(X(7)))-(Dist(2)*K(2)*X(2))/(2*sqrt(X(7))) Dist(2)*K2()*sqrt(X(7))-K(1) 0 0 0;
    0 0 -(K(6)*U(3))/(2*sqrt(X(3))) 0 0 0 0 0 0 0 0;
    K(15)*U(1)/X(3) 0 (X(4)*(K(7)+K(15)*U(1))-K(7)*X(6)-K(15)*U(1)*X(1))/(X(3)^2) -K(7)/X(3)-K(15)*U(1)/X(3) 0 K(7)/X(3) 0 0 0 0 0;
    0 0 0 0 -K(3)*U(2) K(3)*U(2) 0 0 0 0 0;
    0 0 0 0 K(3)*U(2) -K(5)-K(3)*U(2) 0 0 0 0 0;
    0 0 K(8)*U(3)/(2*sqrt(X(3))) 0 0 0 -Dist(2)*K(9)/(2*sqrt(X(7))) 0 0 0 0;
    0 0 (K(8)*U(3)*(X(4)-X(8)))/(2*sqrt(X(3))*X(7)) K(8)*U(3)*sqrt(X(3))/X(7) 0 0 K(8)*U(3)*sqrt(X(3))*(X(8)-X(4))/(X(7)^2) -K(8)*U(3)*sqrt(X(3))/X(7) 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 -K(11)*U(4) K(11)*U(4);
    0 0 0 0 0 0 0 0 K(12) K(11)*U(4) -K(12)-K11)*U(4)];

b = [Dist(1)*K(12)-K(12)*X(1) 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    K(15) 0 -K(6)*sqrt(X(3)) 0 0 0 0 0;
    K(15)*(X(1)-X(4))/X(3) 0 0 0 0 0 0 0;
    0 K(3)*X(6)-K(3)*X(5) 0 0 K(4) 0 0 0;
    0 K(3)*X(5)-K(3)*X(6) 0 K(5) 0 0 0 0;
    0 0 K(8)*sqrt(X(3)) 0 0 0 0 0;
    0 0 K(8)*sqrt(X(3))*(X(4)-X(8))/X(7) 0 0 0 K(10) 0;
    0 0 0 0 0 K(13) -K(13) K(13);
    0 0 0 K(11)*(X(11)-X(10)) 0 K(4) 0 0;
    0 0 0 K(11)*(X(10)-X(11)) 0 0 0 0];

b = [b dxdt];

% Make state indices consistent with x1 = h, x2 = T_T
% ix = [2,1];
% 
% a = a(ix,ix);
% b = b(ix,:);

function [A, B] = getDiscrete(a, b, Ts)
% Convert to discrete time
A = expm(a*Ts);
nx = size(b,1);
n = 4;  % Number of points for Simpson's Rule, an even integer >= 2.
% Use Simpson's rule to compute integral(0,Ts){expm(a*s)*ds*b}
h = Ts/n;
Ai = eye(nx) + A;        % First and last terms;
Coef = 2;
for i = 1:n-1
    if Coef == 2
        Coef = 4;
    else
        Coef = 2;
    end
    Ai = Ai + Coef*expm(a*i*h);     % Intermediate terms
end
B = (h/3)*Ai*b;
        