clc; clear all; close all;
% Lars Reiter Nielsen 
% Pertubation og stabilitet af 5 node topologi
% model: Runge kutta

% Grundparametre
T = 400; % tid for simuleringen (sekunder)
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

% Koblingsstyrke og forbindelsesmatrix
pmax = 2.13; % antager samme styrke 
A = [[0 1 0 1 0]; [1 0 1 0 0]; [0 1 0 1 0]; [1 0 1 0 1]; [0 0 0 1 0]];

 
% vægte i RK4
q = [1/6; 1/3; 1/3; 1/6];

% faseforskels vektor til plotning
Dtheta_12 = [];

% define phase order parameter list and pmaxes list
Rinf = [];
times = [];
pmaxes = [];

Ps = P; % set power list

for p=2.09:0.0025:2.16
    % initialiser
    R = [];
    pmaxes(end+1) = p;
    F = [0; 0; 0; 0; 0;];
    W = [0; 0; 0; 0; 0;];
for t = 0:dt:T

% opdater approksimationer til f(t, theta, omega) og g(t, theta, omega)
k1 = Ps - alpha*W + oscil(p*A,F);
l1 = W;
k2 = Ps - alpha*(W+k1*dt/2)+oscil(p*A,F+l1*dt/2);
l2 = W+k1*dt/2;
k3 = Ps - alpha*(W+k2*dt/2)+oscil(p*A,F+l2*dt/2);
l3 = W+k2*dt/2;
k4 = Ps - alpha*(W+k3*dt)+oscil(p*A,F+l3*dt);
l4 = W+k3*dt;

K = [k1 k2 k3 k4]; % Nx4 matrix med k_1(i) = f(t, w_i, theta_i) etc.
L = [l1 l2 l3 l4];

F = F + dt*(L*q); % opdater theta (fasen), hvor F(1) er theta_1
W = W + dt*(K*q); % opdater w (fasehastigheden), W(1) er w_1

    if t > T/2
        R(end + 1) = 1/N*(sum(exp(1i*F)));
        times(end+1) = t;
    end
end
Rinf(end+1) = real(sum(R))/length(R);

% recognize critical value
if p > 2.09 && (Rinf(end)-Rinf(end-1)) > 0.1
    x_p = [(p+(p-0.0025))/2, (p+(p-0.0025))/2];
    y_p = [-0.01,1];
end
end

hold on
plot(x_p,y_p, 'r-','LineWidth', 2)
plot(pmaxes, Rinf, 'b--o')
title('Critical coupling strength')
xlabel('K (coupling strength)')
ylabel('r_{\infty}          ', 'rot',0)
axis([2.088 2.162 -0.01 1])
set(gca, 'XTick',[2.09 2.10 2.11 x_p(1) 2.14 2.15 2.16]);
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

