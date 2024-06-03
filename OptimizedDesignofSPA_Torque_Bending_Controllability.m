%% Optimization of Soft Actuator
% Date: 06/07/2023
clear all;
close all;
clc;
% This m file aims to find the optimal dimensional parameters of soft
% pneumatic actuators(SPA). Two kinematic models are obtained by mechanical
% analysis - "Pressure-to-Toruqe model" and "Pressure-to-Bending model".
% Both kinematic models are served as the objective function in optimization 
% algorithm to enhance both torque and bending of SPA.
% A dynamical equation is built to served as a constraint which places the 
% "natural frequency" in a desire range.



%% Define variables of design space
% a is the distance between neutral surface and bottom
a = optimvar('a', 'Lowerbound', 0.002, 'UpperBound', 0.004);
% b is the distance between neutral surface and top
b = optimvar('b', 'LowerBound', 0.014, 'UpperBound', 0.024);
% t is the wall thickness of the chamber room in cross-sectional area
t = optimvar('t', 'LowerBound', 0.0015, 'UpperBound', 0.003);
% w is the width of the cross-sectional of soft actuator
w = optimvar('w', 'LowerBound', 0.010, 'UpperBound', 0.030);
% P is the input pressure and here is set as a constant
P = optimvar('P');
% wn is the upper bound of natural frequency
wu = 3.5;
% wl is the lower bound of natural frequency
wl = 2.5;


%% Create objective function
% Pressure-to-torque model
f1 = 0.5*P*(b-t)*(b-t)*(w-2*t);
f2 = ((b-t)*(w-2*t)*P)/(a*w+w*t+2*b*t-2*t^2);
f3 = 0.5*(w*a^2-(w*t)*(2*b-t)-2*t*(b-t)*(b-t));
C = (f1+f2*f3)/P; % This parameter will be used in Pressure-to-Bending model
torque = (f1+f2*f3);
% If we set torque as "0", the algo only optimizes the bending angle

% Pressure-to-bending model 
I = (1/2)^(1+n)*(1/(2+n))*w*(b+a)^(2+n); % moment of inertia
A = (w-2*t)*(b-t);  % cross-sectional area of chamber only
Aw = (w*(a+b))-A;   % cross-sectional area of SPA minus chamber area
L = 0.094;    % Lenght of SPA
n = 1.0;      % n varies with material properties; n=1 for linear model; n>1 for nonlinear model
E = 340000;   % Young's modulus of DragonSkin 20
theta = (((n/(n+1))^n*(C*P/(E*I))*(L+(L*P*A/(Aw*E)))^n))^(1/n)*180/pi;
% If we set theta as "0", the algo only optimizes the torque

%% Determine the optimization function
% objective function
prob = optimproblem("Objective", -(torque/0.36 * theta/250));
options = optimoptions('fmincon', 'Algorithm','interior-point','Display','iter','ConstraintTolerance',1e-12);
% options = optimoptions('fmincon', 'Algorithm','sqp','Display','iter','ConstraintTolerance',1e-12);
% The torque and bending angle are in different units.
% Thus, torque is normalized by 0.36 N-m and theta is normalized by 250 deg
% to balance both performance metrics

% Define additional constraints
prob.Constraints.cons1 = P == 0.15*1000000;
prob.Constraints.cons2 = a + b <= 0.024;  % a + b is the height of cross-sectional area
prob.Constraints.cons3 = a + b >= 0.018;
% Constrains to place the natural frequency
prob.Constraints.cons4 = w*(a+b)^(n+2) >= ((n/(n+1))^n*2^(n+1)*(n+2)*0.36*wl^2*L^(n+1))/E; % Wn >= 2.5
prob.Constraints.cons5 = w*(a+b)^(n+2) <= ((n/(n+1))^n*2^(n+1)*(n+2)*0.36*wu^2*L^(n+1))/E; % Wn <= 3.5

% Define the initial conditions for the algorithms to search
x0.a = 0.003;
x0.b = 0.018;
x0.t = 0.0025;
x0.w = 0.025;
x0.P = 0.15*1000000;
% Select the solver as "interior-point"
opts = optimoptions(@fmincon,'Algorithm','interior-point')
[sol, fval, exitflag, output] = solve(prob, x0, 'Options', opts);
sol
optimzed_val = -fval

%% To check if the natural frequency is in desired range
n = 1.0; % n dependes on the soft materials (1.0~2.0)
naturalFreq = sol.w*(sol.a + sol.b)^(n+2)
