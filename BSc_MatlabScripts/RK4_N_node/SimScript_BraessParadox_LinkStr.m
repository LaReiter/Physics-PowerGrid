clc; clear all; close all;
%%% Simulation of the power grid kuramato model %%%
%%% Checking for Braess Paradox by STRENGTHENING present links in 30 node
%%% topologies with a danish network structure

% ground parameters
dt = 1e-1; % timestep
alpha = 0.1; % dissipation parameter
totalTime = 400; % time duration
N0 = 30; % number of oscillators
pmax = 3; % standard line capacity


% (1) Load initial conditions
init = load('TopMat_InitVal.txt');
M = load('GridConMat.txt');

% Investigate Braess Paradox on random topology based on danish maps
% structure.
% Outer loop: Repeat procedure: create new topology
% First (while)loop: Define random topology matrix and assert stability.
% Second (for) loop: run over all possible links prone to BP

Topologies = []; % contains number of nodes in the topology for each iteration
Braess = []; % contains fraction of links showing BP for each topology

for j = 1:40

Rinf = 0;

while Rinf < 0.1
    % (1) Load initial conditions
    edges = randi(15);
    A = CreateAdj2(N0,pmax,2*N0-edges,M);
    
    init = load('initval.txt');
    
    W0 = transpose(init(1,:)); % initial condition near fixpoint
    T0 = transpose(init(2,:)); % initial condition near fixpoint
    P0 = init(3,:); % power configuration
    
    s1 = abs(sum(P0));
    
    P0([11,12,13,22]) = [3/20*s1 3/20*s1 3/20*s1 3/20*s1];
    P0([2,6,10,19,28]) = [2/25*s1 2/25*s1 2/25*s1 2/25*s1 2/25*s1];
    P0 = transpose(P0);
    
    [times,timesR, theta, power, R] = NnodeSim(totalTime, dt, alpha, A, P0, W0, T0);

    Rinf = sum(R)/length(R); 
end

[Kc, Rinf]=findKc(totalTime,dt,alpha,A,P0,W0,T0); % find critical Coupling Strength

idx = find(triu(A)); % find non-zero elements in upper A (potential new links, prone to BP)
links = length(idx); % total number of links in topology prone to BP
BraessCounter = 0;

% investigate BP on each link 
for i = transpose(idx)
    
    A0 = A;
    A0(i) = A0(i)+2; % add capacity
    A0 = triu(A0)+triu(A0).'; % symmetric matrix
    
   [Kc_BP, Rinf_BP] = findKc(totalTime, dt, alpha, A0, P0, W0, T0);
    
    if (Kc_BP > Kc) && (Rinf_BP > 0.1);
        BraessCounter =  BraessCounter + 1;
    end

end

Braess(end+1) = BraessCounter/links;
Topologies(end+1) = j;

disp(j); % keep track of time
end
