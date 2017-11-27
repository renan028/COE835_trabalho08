%----------------------------------------------------------------------
%
%  COE-835  Controle adaptativo
%
%  Script para simular o trabalho 7
%
%  Backstepping  :  n  = 2     Second and third order plant
%                   n* = 2     Relative degree
%                   np = 3     Adaptive parameters
% Caso com observador de ordem reduzida
%----------------------------------------------------------------------

function dx = backstepping_red(t,x)

global A B thetas N c1 c2 d1 d2 Gamma gamma kp a w e1;

y           = x(1:2);
theta       = x(3:5);
lambda      = x(6);
eta         = x(7);
rho         = x(8);

%% Input
yr = a(1) * sin(w(1)*t) + a(2) * sin(w(2)*t);
dyr = a(1) * w(1) * cos(w(1)*t) + a(2) * w(2) * cos(w(2)*t);
ddyr = -a(1) * w(1)^2 * sin(w(1)*t) - a(2) * w(2)^2 * sin(w(2)*t);

Phi = [-y(1) 0;0 -y(1)];

%% Variables 1
xi = -N^2 * eta;
Xi = -[N*eta eta];
v1 = lambda(1);
omega_bar = [0, (Xi - y(1)*e1')]';
omega = [v1, (Xi - y(1)*e1')]';

%% Z
z1 = y(1) - yr;
alpha_bar = -c1*z1 - d1*z1 - xi - omega_bar'*theta;
alpha_1 = rho * alpha_bar;
z2 = v1 - rho*dyr - alpha_1;

%% Filtro eta
deta = N*eta + y(1);

%% dalpha/dt
dady = rho * (- c1 - d1 + [0,e1']*theta); 
dadeta_deta = rho * (N^2 * deta + [0, N*deta, deta]*theta);
dadyr = rho*(c1 + d1);
dadtheta = - rho * omega_bar';
dadrho = -(c1 + d1)*(y(1) - yr) - xi - omega_bar'*theta;


%% Variables 2
tau_1 = (omega - rho*(dyr + alpha_bar)*[e1',0]')*z1;
tau_2 = tau_1 - z2 * (dady * omega); 

%% Atualização dos parâmetros
dtheta = Gamma * tau_2;
drho = - gamma * z1 * sign(kp) * (dyr + alpha_bar);
beta = N*v1 + dady * (xi + omega'*theta) + ...
    dadeta_deta + dadyr * dyr + (dyr + dadrho) * drho;
u = -c2*z2 + beta + rho*ddyr + dadtheta*dtheta - d2*z2*(dady)^2 - ...
    theta(1)*z1;

%% Filtros
dlambda = N*lambda + u;


%% Planta
F = [B*u Phi];
dy = A*y + F*thetas;

%% Translation
dx = [dy' dtheta' dlambda' deta' drho]';    
