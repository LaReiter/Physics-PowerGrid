clc; clear all; close all;
% Lars Reiter Nielsen 
% 5 node Braess paradox. Parametre: alpha = 1. P = [1,-2,-2,2,1].Pmax =2.13
% styrk kapacitet på forbindelse 4-5 for at opnå braess paradox
% Topologi med 4 generatorer og 4 forbruger
% model: Runge kutta

% Grundparametre
T = 100; % tid for simuleringen (sekunder)
dt = 1e-1; % BØR VÆLGES SMART
Omega = 50; % standard frekvens (hz), 50-60hz
alpha = 1; % dissipations parameter
tracktime = 0; % timedelay on data (time/phase difference) collection
N = 5; % antal oscillatorer


% Effekt for G_1 (P1), M (P2) og G_2 (P3)
P1 = 1; P2 =-2 ; P3 = -2; P4 = 2; P5 = 1;
P = [P1; P2; P3; P4; P5;];

% frekvenser: theta_i er fasen af oscillator i og w_i er fasehastigheden

% % theta_1 = Omega + randsample([-1,1],1)*rand(1)*scale; % afvigelse fra standard frekvensen med +-(0..2)
% % theta_2 = Omega + randsample([-1,1],1)*rand(1)*scale; % randsample([-1,1],1) eksluderer muligheden 0
% % theta_3 = Omega + randsample([-1,1],1)*rand(1)*scale;
F = [0; 0; 0; 0; 0;];

W = [0; 0; 0; 0; 0;];

% Koblingsstyrke og forbindelsesmatrix
pmax = 2.13; % antager samme styrke 
A = pmax*CreateAdj(N);

% vægte i RK4
q = [1/6; 1/3; 1/3; 1/6];

% faseforskels vektor til plotning
Dtheta_12 = [];

% theta cell array og tids cell array til multiplot
Theta = {};
times = [];

for t = 0:dt:T
Ps = P;

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
title('Normal operation')
xlabel('time (a.u.)')
ylabel('\phi(t)    ','rot',0)
legend('\phi_1', '\phi_2', '\phi_3', '\phi_4', '\phi_5')
hold off


