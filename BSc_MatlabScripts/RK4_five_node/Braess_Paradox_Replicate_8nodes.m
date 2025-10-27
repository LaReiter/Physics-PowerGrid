clc; clear all; close all;
% Lars Reiter Nielsen 
% REPLIKATION AF BRAESS PARADOX (Se Braess Paradox artikel, s.4)
% Topologi med 4 generatorer og 4 forbruger
% model: Runge kutta

% Grundparametre
T = 40; % tid for simuleringen (sekunder)
dt = 1e-1; % BØR VÆLGES SMART
Omega = 50; % standard frekvens (hz), 50-60hz
alpha = 1; % dissipations parameter
tracktime = 0; % timedelay on data (time/phase difference) collection
N = 3; % antal oscillatorer


% Effekt for G_1 (P1), M (P2) og G_2 (P3)
P1 = -1; P2 = 1; P3 = 1; P4 = 1; P5 = -1; P6 = -1; P7 = 1; P8 = -1;
P = [P1; P2; P3; P4; P5; P6; P7; P8];

% frekvenser: theta_i er fasen af oscillator i og w_i er fasehastigheden

% % theta_1 = Omega + randsample([-1,1],1)*rand(1)*scale; % afvigelse fra standard frekvensen med +-(0..2)
% % theta_2 = Omega + randsample([-1,1],1)*rand(1)*scale; % randsample([-1,1],1) eksluderer muligheden 0
% % theta_3 = Omega + randsample([-1,1],1)*rand(1)*scale;
F = [0; 0; 0; 0; 0; 0; 0; 0;];

W = [0; 0; 0; 0; 0; 0; 0; 0;];

% Koblingsstyrke og forbindelsesmatrix
pmax = 1.03; % antager samme styrke
A = zeros(length(F), length(F)); % 1 betyder kobling mellem oscillatorpar (i,j)
A(1,2) = 1; A(2,1) = 1; A(1,7) = 1; A(7,1) = 1; A(1,6) = 1; A(6,1) =1;
A(2,3) = 1; A(3,2) =1; A(3,7) = 1; A(7,3) =1; A(3,4) = 1; A(4,3) = 1; 
A(4,8) =1; A(8,4) =1; A(4,5) =1; A(5,4) =1; A(5,6) =1; A(6,5) =1; A(6,8) =1;
A(8,6) =1;
A = pmax*A;
A(4,3) =10; A(3,4) =10;

 
% vægte i RK4
q = [1/6; 1/3; 1/3; 1/6];

% faseforskels vektor til plotning
Dtheta_12 = [];

% theta cell array og tids cell array til multiplot
Theta = {};
times = [];

% Effekt pertubation (forøg forbruget)
% dP = P + [-4.9; 0; 0];   % NB kritisk punkt -4.9 dP

for t = 0:dt:T
Ps = P;
% Pertubation
% if (t < 42) && (40 < t)
%     Ps = dP;
% else
%     Ps = P;
% end

% opdater approksimationer til f(t, theta, omega) og g(t, theta, omega)
k1 = Ps - alpha*W + oscil(A,F);
l1 = W;
k2 = Ps - alpha*(W+k1*dt/2)+oscil(A,F+l1*dt/2);
l2 = W+k1*dt/2;
k3 = Ps - alpha*(W+k2*dt/2)+oscil(A,F+l2*dt/2);
l3 = W+k2*dt/2;
k4 = Ps - alpha*(W+k3*dt)+oscil(A,F+l3*dt);
l4 = W+k3*dt;

K = [k1 k2 k3 k4]; % Nx4 matrix med k_1(i) = f(t, w_i, theta_i) etc.
L = [l1 l2 l3 l4];

F = F + dt*(L*q); % opdater theta (fasen), hvor F(1) er theta_1
W = W + dt*(K*q); % opdater w (fasehastigheden), W(1) er w_1

% theta difference til plot
Dtheta = [F - F(1) F - F(2) F - F(3)];

% theta til plot
Theta{end+1} = F;

% tid til plot
times(end+1) = t;
end

hold on
plot(times, cell2mat(Theta))
hold off


% % hold on
% % [hAx, hLine1, hLine2] = plotyy(times, Dtheta_12, times, pmax*sin(Dtheta_12));
% % title('Pertubation behaviour');
% % xlabel('time (a.u.)');
% % ylabel(hAx(1), '\theta_1-\theta_2');
% % ylabel(hAx(2), 'Supplied power');
% % set(hLine1, 'LineStyle', '--');
% % set(hLine2, 'LineStyle', '-');
% % set(hLine1, 'LineWidth', 3);
% % hold off

