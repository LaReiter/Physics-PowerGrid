clc; clear; close all;
format compact; % less spacing

% Test one node basin stability of arbitrary network

%% Initialization
% N = 20; % number of oscillators
split = 0.5; % fraction of consumers
uniform = true; % make distribution and generation equal for both consumers and generators
Pmax = 5; % maximum consumer/generator power demand for power vector
d = 2.7; % average degree
Time = 100; % runtime
dt = 0.1; % timestep
alpha = 0.1; % dissipation parameter
Trials = 50; % number of pertubation trials in onenodeBasinStability
% Omega0 = ones(N,1); % initial vector for phase velocity
% Theta0 = ones(N,1); % initial vector for phase
%%

% A = createAdjMat(N,d);
% % CHECKPOINT
% disp('checkpoint 1');
% 
% P = createPowerVec(A,split,Pmax,uniform);
% % CHECKPOINT
% disp('checkpoint 2');
[A,P] = createTopology(4,'R');
N = length(A);
Omega0 = ones(N,1); % initial vector for phase velocity
Theta0 = ones(N,1); % initial vector for phase

K = round(max(abs(P)))*1.5;  % critical coupling strength (strong coupling!)
% CHECKPOINT
disp('checkpoint 3');

[bool, X] = simTraj(2,Time,dt,alpha,K*A,P',Omega0,Theta0);
% CHECKPOINT
disp('checkpoint 4');

if bool == 1
    onenodeBasinStability(K*A, X(:,2),X(:,1),P',alpha,dt,Time,Trials,'Uniform')
end