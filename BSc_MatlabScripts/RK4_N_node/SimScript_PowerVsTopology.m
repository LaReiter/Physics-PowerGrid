clc; clear all; close all;
%%% Simulation of the power grid kuramato model %%%

% ground parameters
dt = 1e-1; % timestep
alpha = 0.1; % dissipation parameter
totalTime = 400; % time duration
N0 = 30; % number of oscillators
P_its = 8; % max power iteration
Rep_its = 15; % number of repetitions of each topology


% % % (1) Load initial conditions
% % A = load('GridConMat.txt');
% % init = load('initval.txt');
% % 
% % W0 = transpose(init(1,:)); % initial condition near fixpoint
% % T0 = transpose(init(2,:)); % initial condition near fixpoint
% % P0 = transpose(init(3,:)); % power configuration

finalRinf = zeros(1,P_its*2); 

for k = 1:1:Rep_its

Rinf = []; % phase order long time average
pmaxes = []; % pmaxes

for p=0.5:0.5:P_its
    
pmaxes(end+1)= p;

%%% (2)  define initial conditions
A = p*CreateAdj(N0);

% define W0 and T0
W0 = zeros(N0,1);
T0 = zeros(N0,1);

% define P0; algorithm to initialize random power vector
P0 = rand(N0,1); % create random power units

g1 = randi(N0); g2 = randi(N0); g3 = randi(N0); g4 = randi(N0); g5 = randi(N0); % pick 5 random generators
s = sum(P0) - P0(g1) - P0(g2) - P0(g3) - P0(g4) - P0(g5); % define the total consumed power s

w1 = rand(1); w2 = (1-w1)*rand(1); w3 = (1-(w1+w2))*rand(1); w4 = (1-(w1+w2+w3))*rand(1); % define weights for generator power distribution
w5 = 1 - (w1+w2+w3+w4);

P0 = -P0; % ensure that consumers drain power
P0(g1) = s*w1; P0(g2) = s*w2; P0(g3) = s*w3; P0(g4) = s*w4; P0(g5) = s*w5; % make the random generators supply power with given weights



[Rinf] = NnodeSim_CritVal(totalTime, dt, alpha, 3*A, P0, W0, T0);

% hold on
% plot(times, cell2mat(theta))
% plot(timesR, R)

disp(p); % tracker

end

finalRinf = finalRinf + Rinf;

end

hold on
title('Stability uniformity criteria')
xlabel('p_{max}','FontWeight' ,'bold');
ylabel({'r_{\infty}', '(20-sim average)'},'FontWeight', 'bold', 'Color', 'b');
plot(pmaxes, (1/Rep_its)*finalRinf,'bs');
hold off





















%%% UNCOMMENT BELOW CODE %%%

% % % % TWO NODE SIM
% % % A=[[0 1]; [1 0]];
% % % pmax = 3; % maximum coupling strength
% % % W0 = [-0.0001;0.0001]; % initial condition near fixpoint
% % % T0 = [-29.3596;-29.0202]; % initial condition near fixpoint
% % % P0 = [-1;1]; % power configuration
% % % dt = 1e-1; % timestep
% % % alpha = 0.1; % dissipation parameter
% % % totalTime = 200; % time duration
% % % dP = [-2.5; 0]; % power pertubation
% % % idx = [1,2]; % relevant oscillator connection
% % % 
% % % [times, theta, power] = NnodeSim(totalTime, dt, alpha,dP, pmax*A, P0, W0, T0,idx);
% % % 
% % % hold on
% % % plot(times,theta)
% % % plot(times, power)
% % % hold off


% % % % THREE NODE SIM
% % % A = [[0 1 1]; [1 0 0]; [1 0 0]]; % coupling matrix
% % % pmax = 3; % maximum coupling strength
% % % W0 = [-0.0001;0.0001;0.0001]; % initial condition near fixpoint
% % % T0 = [-29.3596;-29.0202;-29.0202]; % initial condition near fixpoint
% % % P0 = [-2;1;1]; % power configuration
% % % dt = 1e-1; % timestep
% % % alpha = 0.1; % dissipation parameter
% % % totalTime = 200; % time duration
% % % dP = [-2; 0; 0]; % power pertubation
% % % idx = [1,2]; % relevant oscillator connection
% % % 
% % % [times, theta, power] = NnodeSim(totalTime, dt, alpha,dP, pmax*A, P0, W0, T0,idx);
% % % 
% % % hold on
% % % plot(times,theta)
% % % plot(times, power)
% % % hold off