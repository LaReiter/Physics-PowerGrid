clc; clear all; close all;
%%% Simulation of the power grid kuramato model %%%
%%% Checking if there is a connection between critical coupling strength of topology
%%% and number of links in place
%%% NOTICE: Based on danish frame. Assumed constant coupling strength.

% ground parameters
dt = 1e-1; % timestep
alpha = 0.1; % dissipation parameter
totalTime = 200; % time duration
N0 = 30; % number of oscillators
eps = 1e-2;  % critical coupling strength error margin
repMax = 10; % number of pmax samples pr links
edgesMin = 29; % minimum number of edges for connectivity
edgesMax = 100; % max number of edges
maxIts = 25; % Set max search iterations for steady state
M = load('GridConMat_EqCouplStr.txt');


% Outer loop: Repeat procedure:
% Second loop: Run over all links and investigate Braess paradox on each
% First (while)loop: Define random topology matrix and assert stability.

links = []; % list for number of topology links
pmaxAv = []; % list for pmax average

for edges = edgesMin:2:edgesMax
    
    links(end+1) = edges; % save edges
    pmaxSample = []; % sample pmaxes for random topologies with fixed number of links
    rep = 0;
    deadendcounter = 0; % +1 if no stability can be found in current rep

    while (rep < repMax) && (deadendcounter < maxIts);
        
        Rinf = 0; % phase order parameter (long time limit)
        pmin = 0; % lowest coupling 
        pmax = 10; % highest coupling
        p = pmax; % current coupl str
        pdiff = pmax;
        
        steadystate = 0; % check if steady state is found
 
        % create coupling matrix
        A = CreateAdj(N0,edges,M);
       
        W0 = zeros(N0,1); % initial condition near fixpoint
        T0 = zeros(N0,1); % initial condition near fixpoint
        P0 = -ones(1,N0); % power configuration

        % power consumption. let renewable sources have the heighest weight
        delta = ( 1/5*rand(1) ); % disturbance / uncertainty (som vinden blæser)
        sign = datasample([-1,1],1);
        W_renew = 0.7 + sign*delta; % one sign
        W_plant = 0.3 - sign*delta; % other sign
        s1 = 21; % 21 consumers  
        P0([11,12,13,22]) = [1/4*W_renew*s1 1/4*W_renew*s1 1/4*W_renew*s1 1/4*W_renew*s1];
        P0([2,6,10,19,28]) = [1/5*W_plant*s1 1/5*W_plant*s1 1/5*W_plant*s1 1/5*W_plant*s1 1/5*W_plant*s1];
        P0 = transpose(P0);        

        while (pdiff > eps)

            pp = p;  % save previous coupling point
            p = (pmax+pmin)/2; % new coupling point
            pdiff = abs(p - pp); % difference from previous and current coupl str

            A0 = p*A; % new coupling strength to coupling matrix

            Rinf = NnodeSim_CritVal(totalTime, dt, alpha, A0,P0, W0, T0);

            if Rinf > 0.9 % steady state reached
                pmax = p;
                steadystate = 1;
            else
                pmin = p;
            end

        end
        
        if steadystate
            pmaxSample(end+1) = p; % save sample pmax
            rep = rep + 1;
            deadendcounter = 0;
        else
            deadendcounter = deadendcounter + 1;
        end
        
%         disp(rep)
%             
    end
    
    disp(edges);
    pmaxAv(end+1) = sum(pmaxSample)/repMax;
    
end


