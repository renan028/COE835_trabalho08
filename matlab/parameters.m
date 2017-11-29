clear;
clc;
close all;

global A B thetas N c1 c2 d1 d2 Gamma gamma kp a w e1;

PRINT = true;
% PRINT = false;

%Simulation time
tfinal = 60;

% Unit vectors
e1 = [1 0]';

% System matrix
A = [0 1;0 0];

%% First parameters

kp_1 = 5;
Z_1 = [1];
P_1 = [1 2 1];
thetas_1 = [kp_1 P_1(2) P_1(3)]';

N_1 = -1;

a_1 = [1 1];
w_1 = [1 3];

%Initial conditions
y0_1  = [0 0]';
theta0_1 = [1 0 0]';
eta0_1 = [0]';
lambda0_1 = [0]';
rho0_1 = 1;

%Adaptation gain
Gamma_1 = 10*eye(3);
gamma_1 = 10;
c1_1 = 1;
c2_1 = 1;
d1_1 = 1;
d2_1 = 1;

%% Second parameters

kp_2 = 5;
Z_2 = [1];
P_2 = [1 2 1];
thetas_2 = [kp_2 P_2(2) P_2(3)]';

N_2 = -1;

a_2 = [1 1];
w_2 = [1 3];

%Initial conditions
y0_2  = [0 0]';
theta0_2 = [1 1 1]';
eta0_2 = [0]';
lambda0_2 = [0]';
rho0_2 = 1;

%Adaptation gain
Gamma_2 = 10*eye(3);
gamma_2 = 10;
c1_2 = 1;
c2_2 = 1;
d1_2 = 1;
d2_2 = 1;