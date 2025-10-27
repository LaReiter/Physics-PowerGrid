 function [ X ] = findSteadyState( Time, dt, alpha, A, P, Omega_start, Theta_start)
% Incrementaly checks for steady state. Return X = [Omega, Theta] if found.
% Time:     Period of simulation
% dt:   Timestep
% alpha:    Disipation parameter
% A:    Connection matrix WITH coupling strength
% P:   Power vector. P(i) > 0 generator. P(i) < 0 consumer.
% Omega_start:  Initial guess for Omega vector
% Theta_start:  Initial guess for Theta vector
% Code by Lars Reiter Nielsen @ Student at Copenhagen University

% Initialization
N = numel(Omega_start); % number of oscillators

sampleTime = max(10,round(Time/(20*dt))); % duration where phase order and mean-omega is calculated. NOTE: Should always be longer than the oscillation time.

SteadyState = false; % initially, set steady state to false
thetaOffset = 0; % start by checking if Theta_start is a phase-lock state

eps = 1e-1; % epsilon for 'closeness'

% define F(X)
F = @(X_in, P_in) [(P_in - alpha*X_in(:,1) + sum(A.*sin(transpose(repmat(X_in(:,2),1,N)) - repmat(X_in(:,2),1,N)),2)) X_in(:,1)];

while ~SteadyState && thetaOffset < 51
    
    %     %%% TESTING
%     TimeList = zeros(1,Time/dt);
%     OmegaList = zeros(N,Time/dt);
        
    currentIteration = 1; % keep track of loop iteration
    
    OmegaInfVec = ones(1,sampleTime);  % vector containing <omega(t_i)> at each discrete timestep t_i
    PhaseOrderVec = ones(1,sampleTime); % vector containing phase order parameter r(t_i) at each discrete timestep t_i

    % Try new initial state
    T = Theta_start*thetaOffset;
    W = Omega_start;

    % Define state vector X = (W,Theta)^T 
    X = [W T];

%     counter2 = 1; %TESTING

    for t = 0:dt:Time

        % opdater approksimationer til f(X)
        K1 = F(X, P);
        K2 = F(X+dt/2*K1, P);
        K3 = F(X+dt/2*K2, P);
        K4 = F(X+dt*K3, P);

        % opdater X = (F, W)^T
        X = X + dt/6*(K1+2*K2+2*K3+K4);
        
%         %%% TESTING
%         OmegaList(:,counter2) = X(:,1);
%         TimeList(counter2) = t;
%         counter2 = counter2 +1;
        
        OmegaInfVec(currentIteration) = 1/N*sum(abs(X(:,1)));
        PhaseOrderVec(currentIteration) = 1/N*sum(exp(1i*X(:,2)));

        currentIteration = currentIteration + 1;
        
        if mod(currentIteration,round(sampleTime/dt)) == 0
            currentIteration = 1;
        end
        
    end
    
    if 1/length(OmegaInfVec)*sum(OmegaInfVec) < eps && abs(1/length(PhaseOrderVec)*abs(sum(PhaseOrderVec)) - 1) < 0.2
        SteadyState = true;
        break
    end
    
    % advance dtheta
    thetaOffset = thetaOffset + 5;    

end

if ~SteadyState
    X = 0;
end

% plot(TimeList, OmegaList)

end

