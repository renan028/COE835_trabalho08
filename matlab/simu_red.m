%----------------------------------------------------------------------
%
%  COE-835  Controle adaptativo
%
%  Script para simular o trabalho 8 
%
%  Backstepping  :  n  = 2     Second and third order plant
%                   n* = 2     Relative degree
%                   np = 3     Adaptive parameters
% Com observador de ordem reduzida
%----------------------------------------------------------------------

global A B thetas N c1 c2 d1 d2 Gamma gamma kp a w e1 k;

sim_str = strcat('sim0','_');

%% Simulation 1 (default)
kp = kp_1;
Z = Z_1;
P = P_1;

thetas = thetas_1;

N = N_1;

a = a_1;
w = w_1;
c1 = c1_1;
c2 = c2_1;
d1 = d1_1;
d2 = d2_1;
Gamma = Gamma_1;
gamma = gamma_1;

% Initialization
y0  = y0_1;
theta0 = theta0_1;
eta0 = eta0_1;
lambda0 = lambda0_1;
rho0 = rho0_1;
init = [y0' theta0' lambda0' eta0' rho0]';

[T_1,X_1] = ode23s('backstepping_red',tfinal,init,'');
y_1      = X_1(:,1);
theta_1 =  X_1(:,3:5);
tiltheta_1 = theta_1 - thetas';
modtt_1 = sqrt(sum(theta_1.^2,2));
r_1 = 0;
for i=1:length(a)
    r_1 = r_1 + a(i)*sin(w(i).*T_1);
end
e0_1 =  y_1 - r_1;


% %% Simulation 2 (gamma)
% changed = 1;
% 
% kp = kp_1;
% Z = Z_1;
% P = P_1;
% 
% thetas = thetas_1;
% 
% 
% N = N_1;
% 
% a = a_1;
% w = w_1;
% c1 = c1_1;
% c2 = c2_1;
% d1 = d1_1;
% d2 = d2_1;
% Gamma = Gamma_2;
% gamma = gamma_2;
% 
% % Initialization
% y0  = y0_1;
% theta0 = theta0_1;
% eta0 = eta0_1;
% lambda0 = lambda0_1;
% rho0 = rho0_1;
% init = [y0' theta0' lambda0' eta0' rho0]';
% 
% [T_2,X_2] = ode23s('backstepping_red',tfinal,init,'');
% y_2      = X_2(:,1);
% theta_2 =  X_2(:,3:5);
% tiltheta_2 = theta_2 - thetas';
% modtt_2 = sqrt(sum(theta_2.^2,2));
% 
% r_2 = 0;
% for i=1:length(a)
%     r_2 = r_2 + a(i)*sin(w(i).*T_2);
% end
% e0_2 =  y_2 - r_2;
% run plot_bkst.m;

% %% Simulation 3 (planta)
% changed = 2;
% 
% kp = kp_2;
% Z = Z_2;
% P = P_2;
% 
% thetas = thetas_2;
% 
% 
% N = N_1;
% 
% a = a_1;
% w = w_1;
% c1 = c1_1;
% c2 = c2_1;
% d1 = d1_1;
% d2 = d2_1;
% Gamma = Gamma_1;
% gamma = gamma_1;
% 
% % Initialization
% y0  = y0_1;
% theta0 = theta0_1;
% eta0 = eta0_1;
% lambda0 = lambda0_1;
% rho0 = rho0_1;
% init = [y0' theta0' lambda0' eta0' rho0]';
% 
% [T_2,X_2] = ode23s('backstepping_red',tfinal,init,'');
% y_2      = X_2(:,1);
% theta_2 =  X_2(:,3:5);
% tiltheta_2 = theta_2 - thetas';
% modtt_2 = sqrt(sum(theta_2.^2,2));
% 
% r_2 = 0;
% for i=1:length(a)
%     r_2 = r_2 + a(i)*sin(w(i).*T_2);
% end
% e0_2 =  y_2 - r_2;
% run plot_bkst.m;
% 
% %% Simulation 4 (model)
% changed = 3;
% 
% kp = kp_1;
% Z = Z_1;
% P = P_1;
% 
% thetas = thetas_1;
% 
% 
% N = N_1;
% 
% a = a_2;
% w = w_2;
% c1 = c1_1;
% c2 = c2_1;
% d1 = d1_1;
% d2 = d2_1;
% Gamma = Gamma_1;
% gamma = gamma_1;
% 
% % Initialization
% y0  = y0_1;
% theta0 = theta0_1;
% eta0 = eta0_1;
% lambda0 = lambda0_1;
% rho0 = rho0_1;
% init = [y0' theta0' lambda0' eta0' rho0]';
% 
% [T_2,X_2] = ode23s('backstepping_red',tfinal,init,'');
% y_2      = X_2(:,1);
% theta_2 =  X_2(:,3:5);
% tiltheta_2 = theta_2 - thetas';
% modtt_2 = sqrt(sum(theta_2.^2,2));
% 
% r_2 = 0;
% for i=1:length(a)
%     r_2 = r_2 + a(i)*sin(w(i).*T_2);
% end
% e0_2 =  y_2 - r_2;
% run plot_bkst.m;

%% Simulation 5 (y0)
changed = 5;

kp = kp_1;
Z = Z_1;
P = P_1;

thetas = thetas_1;


N = N_2;

a = a_1;
w = w_1;
c1 = c1_1;
c2 = c2_1;
d1 = d1_1;
d2 = d2_1;
Gamma = Gamma_1;
gamma = gamma_1;

% Initialization
y0  = y0_2;
theta0 = theta0_2;
eta0 = eta0_2;
lambda0 = lambda0_2;
rho0 = rho0_2;
init = [y0' theta0' lambda0' eta0' rho0]';

[T_2,X_2] = ode23s('backstepping_red',tfinal,init,'');
y_2      = X_2(:,1);
theta_2 =  X_2(:,3:5);
tiltheta_2 = theta_2 - thetas';
modtt_2 = sqrt(sum(theta_2.^2,2));

r_2 = 0;
for i=1:length(a)
    r_2 = r_2 + a(i)*sin(w(i).*T_2);
end
e0_2 =  y_2 - r_2;
run plot_bkst.m;