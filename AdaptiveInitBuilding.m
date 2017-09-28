%% Adaptive MPC Control of Nonlinear Chemical Reactor Using Successive Linearization
% This example shows how to use an Adaptive MPC controller to control a
% nonlinear continuous stirred tank reactor (CSTR) as it transitions from
% low conversion rate to high conversion rate.
%
% A first principle nonlinear plant model is available and being linearized
% at each control interval. The adaptive MPC controller then updates its
% internal predictive model with the linearized plant model and achieves
% nonlinear control successfully.

% Copyright 1990-2014 The MathWorks, Inc.

%% About the Continuous Stirred Tank Reactor
% A Continuously Stirred Tank Reactor (CSTR) is a common chemical system in
% the process industry. A schematic of the CSTR system is:
%
% <<../mpc_cstr.png>>

%% 
% This is a jacketed non-adiabatic tank reactor described extensively in
% Seborg's book, "Process Dynamics and Control", published by Wiley, 2004.
% The vessel is assumed to be perfectly mixed, and a single first-order
% exothermic and irreversible reaction, A --> B, takes place.  The inlet
% stream of reagent A is fed to the tank at a constant volumetric rate.
% The product stream exits continuously at the same volumetric rate and
% liquid density is constant.  Thus the volume of reacting liquid is
% constant.
%
% The inputs of the CSTR model are:
%
% $$ \begin{array} {ll}
% u_1 = CA_i \; & \textnormal{Concentration of A in inlet feed
% stream} [kgmol/m^3] \\
% u_2 = T_i \; & \textnormal{Inlet feed stream temperature} [K] \\
% u_3 = T_c \; & \textnormal{Jacket coolant temperature} [K] \\
% \end{array} $$
%
% and the outputs (y(t)), which are also the states of the model (x(t)), are:
%
% $$ \begin{array} {ll}
% y_1 = x_1 = CA \; & \textnormal{Concentration of A in reactor tank} [kgmol/m^3] \\
% y_2 = x_2 = T \; & \textnormal{Reactor temperature} [K] \\
% \end{array} $$
%
% The control objective is to maintain the concentration of reagent A, $CA$
% at its desired setpoint, which changes over time when reactor transitions
% from low conversion rate to high conversion rate.  The coolant
% temperature $T_c$ is the manipulated variable used by the MPC controller
% to track the reference as well as reject the measured disturbance arising
% from the inlet feed stream temperature $T_i$.  The inlet feed stream
% concentration, $CA_i$, is assumed to be constant. The Simulink model
% |mpc_cstr_plant| implements the nonlinear CSTR plant.
%
% We also assume that direct measurements of concentrations are unavailable
% or infrequent, which is the usual case in practice.  Instead, we use a
% "soft sensor" to estimate CA based on temperature measurements and the
% plant model.

%% About Adaptive Model Predictive Control
% It is well known that the CSTR dynamics are strongly nonlinear with
% respect to reactor temperature variations and can be open-loop unstable
% during the transition from one operating condition to another.  A single
% MPC controller designed at a particular operating condition cannot give
% satisfactory control performance over a wide operating range.
%%
% To control the nonlinear CSTR plant with linear MPC control technique,
% you have a few options:
%
% * If a linear plant model cannot be obtained at run time, first you need
% to obtain several linear plant models offline at different operating
% conditions that cover the typical operating range.  Next you can choose
% one of the two approaches to implement MPC control strategy:
%
% (1) Design several MPC controllers offline, one for each plant model.  At
% run time, use Multiple MPC Controller block that switches MPC controllers
% from one to another based on a desired scheduling strategy.  See
% <docid:mpc_examples.example-ex20937898>
% for more details.  Use this approach when the plant models have different
% orders or time delays.
%
% (2) Design one MPC controller offline at the initial operating point.  At
% run time, use Adaptive MPC Controller block (updating predictive model at
% each control interval) together with Linear Parameter Varying (LPV)
% System block (supplying linear plant model with a scheduling strategy).
% See <docid:mpc_examples.example-ex18818571> for more details.  Use
% this approach when all the plant models have the same order and time
% delay.
%
% * If a linear plant model can be obtained at run time, you should use
% Adaptive MPC Controller block to achieve nonlinear control.  There are
% two typical ways to obtain a linear plant model online: 
%
% (1) Use successive linearization as shown in this example.  Use this
% approach when a nonlinear plant model is available and can be linearized
% at run time.
%
% (2) Use online estimation to identify a linear model when loop is closed.
% See <docid:mpc_examples.example-ex24638454> for more details.  Use this
% approach when linear plant model cannot be obtained from either an LPV
% system or successive linearization.

%% Obtain Linear Plant Model at Initial Operating Condition
% To linearize the plant, Simulink(R) and Simulink Control Design(R) are
% required.
if ~mpcchecktoolboxinstalled('simulink')
    disp('Simulink(R) is required to run this example.')
    return
end
if ~mpcchecktoolboxinstalled('slcontrol')
    disp('Simulink Control Design(R) is required to run this example.')
    return
end
%%
% To implement an adaptive MPC controller, first you need to design a MPC
% controller at the initial operating point where CAi is 10 kgmol/m^3, Ti
% and Tc are 298.15 K.
%%
%%%%% Definition of Constants %%%%%
g = 9.81 * 1e2 * 3600^2;           % Gravitational constant [cm h^-2]
K_valve = 0.865;    % Valve coefficient conversion constant [to SI units]
rho_water = 1000 * 1e-6;   % Density of water [kg cm^-3]
rho_oil = 800 * 1e-6;      % Density of heating oil [kg cm^-3]
Cp_water = 4180 / 3600;    % Specific heat capacity of water [Wh kg^-1 K^-1]
Cp_oil = 1700 / 3600;      % Specific heat capacity of heating oil [Wh kg^-1 K^-1]
lambda = sqrt(rho_water*g);
Q_oil_tst = 0.2*3600*1e4;  % Heating oil flowrate from TST
Q_oil_tc = 0.2*3600*1e6;   % Heating oil flowrate to PHT
V_tube_tc = 5e6;           % Heat transfer tube volume 
V_tube_tst = 2e6;          % Heat transfer tube volume
V_tst = 30e6;              % Tank volume
A_dhw = 4e2;               % Tank Area
A_s_tst = 0.001e4;         % Heating tube surface area
U_tst = 1275e-4;           % Overall heat transfer coefficient
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
C_v_dhw = 0.0552;

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

%% Definition of inputs
initial_inputs = [0.126, 0.126, 0.126, 0.126, 0.126, 0.126, 0.126, 0.126];
% Create operating point specification.
plant_mdl = 'buildingModel';
op = operspec(plant_mdl);
%%
% Feed concentration is known at the initial condition.
op.Inputs(1).u = initial_inputs';
op.Inputs(1).Known = logical(initial_inputs)';
op.Inputs(2).u = [0, 0]';
op.Inputs(2).Known = [0, 0]';


%%
% Compute initial condition.
[op_point, op_report] = findop(plant_mdl,op); 
%%
% Obtain nominal values of x, y and u.
x0 = [op_report.States(1).x;op_report.States(2).x];
y0 = [op_report.Outputs(1).y;op_report.Outputs(2).y];
u0 = [op_report.Inputs(1).u;op_report.Inputs(2).u];  
%%
% Obtain linear plant model at the initial condition.
sys = linearize(plant_mdl, op_point); 
%%
% Discretize the plant model because Adaptive MPC controller only accepts a
% discrete-time plant model.
Ts = 1;
plant = c2d(sys,Ts);

%% Design MPC Controller
% You design an MPC at the initial operating condition.  When running in
% the adaptive mode, the plant model is updated at run time.
%%
% Specify signal types used in MPC.
plant.InputGroup.MeasuredDisturbances = [];
plant.InputGroup.ManipulatedVariables = [1 2];
plant.OutputGroup.Measured = [1 2];
plant.OutputGroup.Unmeasured = [];
plant.InputName = {'Ti','Tc'};
plant.OutputName = {'T','CA'};
%%
% Create MPC controller with default prediction and control horizons
%set prediction horizon
p = 10;
%set control horizon
m = 2;

mpcobj = mpc(plant,Ts,p,m);
%%
% Set nominal values in the controller
mpcobj.Model.Nominal = struct('X', x0, 'U', u0, 'Y', y0, 'DX', [0 0]);
%%
% Set scale factors because plant input and output signals have different
% orders of magnitude
Uscale = [1 1];
Yscale = [1 1];
mpcobj.MV(1).ScaleFactor = Uscale(1);
mpcobj.MV(2).ScaleFactor = Uscale(2);
mpcobj.OV(1).ScaleFactor = Yscale(1);
mpcobj.OV(2).ScaleFactor = Yscale(2);

mpcobj.MV(1).Min = 0;
mpcobj.MV(1).Max = 10;
mpcobj.MV(2).Min = 0;
mpcobj.MV(2).Max = 10;

mpcobj.OV(1).Min = 0.5;
mpcobj.OV(1).Max = 2;
mpcobj.OV(2).Min = 10;
mpcobj.OV(2).Max = 70;
%%
% Let reactor temperature T float (81.5i.e. with no setpoint tracking error
% penalty), because the objective is to control reactor concentration CA
% and only one manipulated variable (coolant temperature Tc) is available.
mpcobj.Weights.OV = [1 1]; 


%%
% Due to the physical constraint of coolant jacket, Tc rate of change is
% bounded by degrees per minute.
%mpcobj.MV.RateMin = -2;
%mpcobj.MV.RateMax = 2;