    function [Rinf] = NnodeSim_CritVal( time, dt_in, alpha_in, A, P_in, W_in, theta_in)
% Simulates a N_node topology determined by the connection matrix A
% which contains both the coupling strength and connection information:
% (1 or 0) times the coupling strength. T_in is the total period of the
% simulation, P_in is a consume/generation power vector for oscillator
% i=1...N and W_in is the phase angle velocity vector while theta_in is the
% phase angle vector. Both are usually initially zero - Output is the long time limit
% phase order parameter Rinf of the system.

% Grundparametre
dt = dt_in; % BØR VÆLGES SMART
alpha = alpha_in; % dissipations parameter
N = length(W_in); % antal oscillatorer

% frekvenser: theta_i er fasen af oscillator i og w_i er fasehastigheden
T = theta_in;
W = W_in;

% Definer tilstands vektor X = (W,Theta)^T således at systemet hedder 
% dX/dt = F(X)
X = [W T];
unitV = ones(N,1);

% % definer F(X)
f = @(V, P) [(P - alpha*V(:,1) + (A.*sin(transpose(repmat(V(:,2),1,N)) - repmat(V(:,2),1,N)))*unitV) V(:,1)];

Ps = P_in;
Rinf = 0;
    
% tids og faseforskels vektor til plotning
Time = [];
R = [];

for t = 0:dt:time

    % opdater approksimationer til f(X)
    a = f(X, Ps);
    b = f(X+dt/2*a, Ps);
    c = f(X+dt/2*b, Ps);
    d = f(X+dt*c, Ps);
    
    % opdater X = (F, W)^T
    X = X + dt/6*(a+2*b+2*c+d);

    Time(end+1) = t; % update time array

    if t > time/2
        R(end + 1) = 1/N*real((sum(exp(1i*X(:,2)))));
    end

end
    
Rinf = sum(R)/length(R);
    
end
    

