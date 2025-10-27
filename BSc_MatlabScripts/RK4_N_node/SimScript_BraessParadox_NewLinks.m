clc; clear all; close all;
%%% Simulation of the power grid kuramato model %%%
%%% Checking for Braess Paradox by ADDING new links in 30 node
%%% topologies with a danish network structure

% ground parameters
dt = 1e-1; % timestep
alpha = 0.1; % dissipation parameter
totalTime = 400; % time duration
N0 = 30; % number of oscillators
pmax = 3; % standard line capacity
errMarg = 0.005; % tolerance
sampleSize = 100;



% (1) Load initial conditions
init = load('TopMat_InitVal.txt');
M = load('GridConMat.txt');

% Investigate Braess Paradox on random topology based on danish maps
% structure.
% Outer loop: Repeat procedure: create new topology
% First (while)loop: Define random topology matrix and assert stability.
% Second (for) loop: run over all possible links prone to BP

Topology = []; % contains indices for each random topology
Braess = []; % contains fraction of links showing BP for topology i

for j = 1:20

Rinf = 0;

while Rinf < 0.1
    
    edges = randi(15);
    % (1) Load initial conditions
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

idx = find(~(triu(A)+tril(ones(N0,N0)))); % find zero elements in upper A (potential new links, prone to BP)
samp = datasample(idx,sampleSize); % sample size
links = length(idx); % total number of potential new links in topology prone to BP
BraessCounter = 0;

% add new links to potential connections and investigate BP
for i = transpose(samp)
    
    A0 = A;
    A0(i) = pmax; % add link
    A0 = triu(A0)+triu(A0).'; % symmetric matrix
    
   [Kc_BP, Rinf_BP] = findKc(totalTime, dt, alpha, A0, P0, W0, T0);
    
    if (Kc_BP > Kc) && (Kc_BP > 0); % If coupling strength increases desynch provoked
        BraessCounter =  BraessCounter + 1;
        disp(i)
    end

end

Braess(end+1) = BraessCounter/links;
Topology(end+1) = j;

disp(j); % keep track of time
end

% hold on
% title('Stability independence from generator placement')
% xlabel('Generator links (on average)','FontWeight' ,'bold');
% ylabel('r_{\infty}','FontWeight', 'bold', 'Color', 'b');
% plot(links, Rinf,'bs');
% hold off