function [ S ] = onenodeBasinStability( A, Theta0, Omega0, P, alpha, dt, Time, Trials, Dist )
% onenodeBasinStability calculates the basin stability of each node
% Theta0: Theta vector of initial (steady) state. Theta0(i) phase at i
% Omega0: Omega vector of initial (steady) state. Omega0(i) phase vel. at i
% A : Connection matrix with inherent coupling strength
% P : Power vector P. P(i) > 0 generator. P(i)< 0 consumer.
% alpha: Disipation parameter.
% dt: Timestep
% Time: Simulation period
% Trials: Number of samples/repetetions. Basin stability is 'Stable
% trajectories'/ 'Number of trials'.
% Dist: Choose between 'Uniform' or 'Realistic'. The nodewise-power
% peterubations are sampled from the distribution. NOTE: If 'uniform',
% large pertubations are as likely as smaller ones.
% OUTPUT: S(i) contains the basin stability of node i
% Code by Lars Reiter Nielsen @ Student at Copenhagen University

N = length(A);
S = zeros(1,N);

pd1 = makedist(Dist,'lower',-100,'upper',100);  % lower and upper bound for theta distribution
pd2 = makedist(Dist, 'lower', -2*pi, 'upper', 2*pi); % lower and upper bound for omega distribution

for j=1:Trials
    
    PertTheta = random(pd1,1,1);  % draw phase pertubation
    PertOmega = random(pd2,1,1); % draw phase velocity pertubation
    
    for i=1:N % for each node trigger pertubation and check stability
        
        Theta = Theta0; % reset to steadystate
        W = Omega0; % reset to steadystate
        
        Theta(i) = PertTheta; % random initial deviation based on Dist
        W(i) = PertOmega; % random initial deviation based on Dist
       
        % simTraj returns 1 if trajectories converge to steady state
        S(i) = S(i) + simTraj(2, Time, dt, alpha, A, P, W, Theta);
%         disp(i) % keep track of time
    end
    
    % Keep track of runtime. Display trial number in cmd.
    disp(strcat('Trial ', num2str(j)));  
    
end
    
    % calculate fraction of trajectories that returns to steady state for each
    % node 
    S = 1/Trials*S; 


end

