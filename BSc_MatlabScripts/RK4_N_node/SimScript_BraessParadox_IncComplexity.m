clc; clear all; close all;
%%% Simulation of the power grid kuramato model %%%
%%% Checking for Braess Paradox: creating increasingly articulated
%%% networks by adding new small generators to the topology and detecting
%%% the frequency of Braess paradox

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

Nodes = []; % contains number of nodes in the topology for each iteration
Braess = []; % contains fraction of links showing BP for each topology

for oscil = 20:-2:0

Rinf = 0;

while Rinf < 0.1
    % (1) Load initial conditions
    A = CreateAdj2(N0+oscil,pmax,round((N0+oscil)*1.5),M);
    
    W0 = zeros(N0+oscil,1); % initial condition zero vector
    T0 = zeros(N0+oscil,1);  % initial condition zero vector

    % handle power vector
    P0 = -ones(N0,1); % set consumer configuration
    gens = 1/2*rand(oscil,1);
    
    s1 = 21 - abs(sum(gens));
    
    P0([11,12,13,22]) = [3/20*s1 3/20*s1 3/20*s1 3/20*s1];
    P0([2,6,10,19,28]) = [2/25*s1 2/25*s1 2/25*s1 2/25*s1 2/25*s1];
    P0 = [P0;gens];
    
    [Kc, Rinf]=findKc(totalTime,dt,alpha,A,P0,W0,T0); % find critical Coupling Strength

end

idx = find(triu(A)); % find non-zero elements in upper A (potential new links, prone to BP)
links = length(idx) % total number of links in topology prone to BP
BraessCounter = 0;

% investigate BP on each link (by removing it)
for i = transpose(idx)
    
    A0 = A;
    A0(i) = 0; % remove link
    A0 = triu(A0)+triu(A0).'; % symmetric matrix
    
   [Kc_BP, Rinf_BP] = findKc(totalTime, dt, alpha, A0, P0, W0, T0);
    
    if (Kc_BP < Kc) && (Rinf_BP > 0.1);
        BraessCounter =  BraessCounter + 1
    end

end

Braess(end+1) = BraessCounter/links;
Nodes(end+1) = N0+oscil;

disp(oscil); % keep track of time
end



% % % %%% Simulation of the power grid kuramato model %%%
% % % %%% Detecting the frequency of Braess Paradox on increasingly
% % % %%% articulated networks where the number of oscillators increase
% % % 
% % % % ground parameters
% % % dt = 1e-1; % timestep
% % % alpha = 0.1; % dissipation parameter
% % % totalTime = 400; % time duration
% % % N0 = 30; % number of oscillators
% % % pmax = 8; % constant coupling strength
% % % 
% % % % (1) Load initial conditions
% % % init = load('initval.txt');
% % % M = load('GridConMat_EqCouplStr.txt');     
% % % 
% % % % Simulation loop parameters
% % % TopMax = 4; % maximum number of topologies
% % % GenMax = 10; % number of generators atached
% % % RepMax = 10; % decide sample size for each generator attachment loop
% % % 
% % % % stability analysis during the connection of new generator sources
% % % % Outer loop: Run over different topologies
% % % % While loop: determine stable topology
% % % % Inner loop: Run over different # of generators attached (1-10)
% % % % Inner loop: Place generated at random position and repeat 10 times
% % % % (Inner loop) First (while)loop: Define random topology matrix and assert stability.
% % % 
% % % hold on
% % % for top = 1:TopMax
% % % 
% % %     Rinf = 0;
% % % 
% % %     Critvals = []; % contains long time average phase order parameter for topologies
% % %     Generators = []; % list of new generators
% % % 
% % %     % create stable topology
% % %     while Rinf < 0.9
% % %           
% % %             edges = 40 + randi(15); % random ammount of edges (links), minimum of 40, max 55
% % %             A = CreateAdj(N0,edges,M);
% % %   
% % %             W0 = transpose(init(1,:)); % initial condition near fixpoint
% % %             T0 = transpose(init(2,:)); % initial condition near fixpoint
% % %             P0 = init(3,:); % power configuration
% % % 
% % %             s1 = abs(sum(P0));
% % % 
% % %             P0([11,12,13,22]) = [3/20*s1 3/20*s1 3/20*s1 3/20*s1];
% % %             P0([2,6,10,19,28]) = [2/25*s1 2/25*s1 2/25*s1 2/25*s1 2/25*s1];
% % %             P0 = transpose(P0);
% % % 
% % %             [K_c, Rinf] = findKc(totalTime, dt, alpha, pmax*A, P0, W0, T0);
% % % 
% % %     end
% % %     
% % %     % add generators
% % %     for gen = 1:GenMax
% % % 
% % %         GenTemp = [];
% % %         
% % %         % average over added generator
% % %         for j = 1:RepMax % repeat 10 times
% % % 
% % %             P00 = P0; % define template P0
% % % 
% % %             % add new random generators
% % %             Vupper = zeros(N0,gen); % upper part of new matrix
% % %             % make random generator connections to grid
% % %             for l = 1:gen
% % %             % place generator at random link with coupling strength 3
% % %                 Vupper(randi(N0),l) = 3;
% % %             end
% % % 
% % %             Vlower = [Vupper.', zeros(gen,gen)]; % lower part of new matrix
% % % 
% % %             % define power of each (between 0 and 1/2)
% % %             pgens = 1/2*rand(1,gen);
% % % 
% % %             % downgrade already existing generator power so sum(consumed) equals sum(generated) 
% % %             idx = find(P00>0); % generator indices
% % %             normfactor = (sum(P00(idx)) - sum(pgens))/sum(P00(idx)); % normalize
% % %             P00(idx) = P00(idx)*normfactor;
% % %             P00(end+1:end+gen) =  pgens;
% % % 
% % % 
% % %             % Reshape matrix and intial vectors
% % %             A0 = [[A Vupper]; Vlower]; % define reshaped matrix with new connections
% % %             W00 = [W0; zeros(gen,1)];
% % %             T00 = [T0; zeros(gen,1)];
% % % 
% % %             [K_cgen,R_gen] = findKc(totalTime, dt, alpha, pmax*A0, P00, W00, T00);
% % % 
% % %             GenTemp(end+1) = K_cgen;
% % %         end
% % % 
% % % 
% % % 
% % %     Critvals(end+1) = sum(GenTemp)/length(GenTemp);
% % %     Generators(end+1) = gen;
% % % 
% % %     disp(gen); % keep track of time
% % % 
% % %     end
% % % 
% % % 
% % %     plot(Generators, Critvals,'--o')
% % %     
% % % end
% % % 
% % % hold off