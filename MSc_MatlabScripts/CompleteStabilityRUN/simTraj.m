function [ Output1, Output2 ] = simTraj( outputType,Time, dt, alpha, A, P, Omega0, Theta0)
%simTraj Simulate trajectories of grid network. Returns diagnostic output.
%   OutputType:   Choose between value 1 or value 2
%
%       1:  Returns a cell array containing {Omega, Theta, Time} where
%       Theta/Omega is a N x Time/dt+1 and Theta/Omega(i,j) contains the value of 
%       oscillator i at time t_j. Time consist of all
%       t_0,t_1,...,T_(Time/dt+1). In addition creates Omega-Time aswell as
%       Theta-Time plots.
%
%       2:  Returns a pair [bool, X] where bool is the number 1 if trajectories
%       converge to a steady state or the number 0 otherwise. X is the
%       final state (steady state if bool = 1)
%
%   Time:   Simulation runtime
%   dt:     Simulation timestep
%   alpha:   dissipation parameter
%   A:     connection (adjacency) matrix. Entry magnitude constitutes coupling.
%   P:  Power vector. P(i)> 0 defines generators. P(i) < 0 makes consumers.
%   Omega0:    Initial phase vel. vector. Omega0(i) is p.v. of oscillator i
%   Theta0:     Initial phase vector. Theta0(i) is phase of oscillator i

% Code by Lars Reiter Nielsen @ Student at Copenhagen University

%% Initialization

N = numel(Omega0); % number of oscillators
% define N x 2 state vector X = (W,Theta)^T and N x 1 unit vector
X = [Omega0 Theta0];
% define (RHS of differential equation) vector function F(X)
F = @(X_in, P_in) [(P_in - alpha*X_in(:,1) + sum( A.*sin( transpose(repmat(X_in(:,2),1,N)) - repmat(X_in(:,2),1,N) ),2) ) X_in(:,1)];

%% outPutType

% return array containing state vector X = [Omega, Theta] at all times
if outputType == 1
    
    % arrays for plotting 
    Omega_all = zeros(N, Time/dt+1);
    Theta_all = zeros(N, Time/dt+1);
    Time_all = zeros(1, Time/dt+1);
    idx = 1; % dummy index
    
    % simulate trajectories
    for t = 0:dt:Time
        
        Omega_all(:,idx) = X(:,1); % Omega vector at time t
        Theta_all(:,idx) = X(:,2); % Theta vector at time t
        Time_all(idx) = t; % add current time
        

        % update intermediate and full-step approximations to F(X)
        K1 = F(X, P);
        K2 = F(X+dt/2*K1, P);
        K3 = F(X+dt/2*K2, P);
        K4 = F(X+dt*K3, P);

        % Update state vector X
        X = X + dt/6*(K1+2*K2+2*K3+K4);
        
        idx = idx + 1;

    end
    
    % plot omega vs time
    figure
    hold on
    for i = 1:N
        plot(Time_all,Omega_all(i,:), 'DisplayName', ['Oscillator ', num2str(i)])
    end
    legend(gca, 'show')
    xlabel('t')
    ylabel('\omega(t)')
    set(get(gca,'Ylabel'),'rotation',0,'HorizontalAlignment','Right')
    title('Phase velocity vs Time')
    
    % plot theta vs time
    figure
    hold on
    for i = 1:N
        plot(Time_all,Theta_all(i,:), 'DisplayName', ['Oscillator ', num2str(i)])
    end
    legend(gca, 'show')
    xlabel('t')
    ylabel('\theta(t)')
    set(get(gca,'Ylabel'),'rotation',0,'HorizontalAlignment','Right')
    title('Phase vs Time')
    hold off
    
    % save Omega_all, Theta_all, Time_all in cell-array
    Output1 = {Omega_all, Theta_all, Time_all};
    
    
% check if steady state exists
elseif outputType == 2
    
    sampleTime = max(10,round(Time/(20*dt))); % duration where phase order and mean-omega is calculated. NOTE: Should always be longer than the oscillation time.
    currentIteration = 1; % keep track of loop iteration
    
    OmegaInfVec = ones(1,sampleTime);  % vector containing <omega(t_i)> at each discrete timestep t_i
    PhaseOrderVec = ones(1,sampleTime); % vector containing phase order parameter r(t_i) at each discrete timestep t_i
    
    eps = 1e-1; % epsilon for 'closeness'
    
    % simulate trajectories
    for t = 0:dt:Time
        
        % update intermediate and full-step approximations to F(X)
        K1 = F(X, P);
        K2 = F(X+dt/2*K1, P);
        K3 = F(X+dt/2*K2, P);
        K4 = F(X+dt*K3, P);

        % Update state vector X
        X = X + dt/6*(K1+2*K2+2*K3+K4);
       
        OmegaInfVec(currentIteration) = 1/N*sum(abs(X(:,1)));
        PhaseOrderVec(currentIteration) = 1/N*sum(exp(1i*X(:,2)));
        
        currentIteration = currentIteration + 1;

        if mod(currentIteration,sampleTime) == 0
            currentIteration = 1;
        end

    end
    
    if 1/length(OmegaInfVec)*sum(OmegaInfVec) < eps && abs(1/length(PhaseOrderVec)*abs(sum(PhaseOrderVec)) - 1) < 0.2
        Output1 = 1;
        Output2 = X;
    else
        Output1 = 0;
        Output2 = X;
    end
    
    
    
end

end

