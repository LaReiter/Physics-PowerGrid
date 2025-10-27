 function [ Time, TimeR, Theta, Power, R] = NnodeSim( time, dt_in, alpha_in, A, P_in, W_in, theta_in, idx)
% Simulates a N_node topology determined by the connection matrix A
% which contains both the coupling strength and connection information:
% (1 or 0) times the coupling strength. T_in is the total period of the
% simulation, P_in is a consume/generation power vector for oscillator
% i=1...N and W_in is the phase angle velocity vector while theta_in is the
% phase angle vector. Both are usually initially zero. idx is a tuple (i,j)
% which signals that oscillator (i,j) has the relevant phase angle
% difference for plotting - output Dtheta thus contains theta_diff(i,j) for
% all times. Alpha_in is the dissipation parameter. dP is a pertubation
% on oscillator 1.

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


% tids og faseforskels vektor til plotning
Theta = {};
Time = [];
TimeR =[];
Power = [];
R = [];


for t = 0:dt:time
    
Ps = P_in;

% % Pertubation
% if (t < 42) && (40 < t)
%     Ps = P_in + dP;
% end
    
% opdater approksimationer til f(X)
K1 = f(X, Ps);
K2 = f(X+dt/2*K1, Ps);
K3 = f(X+dt/2*K2, Ps);
K4 = f(X+dt*K3, Ps);

% opdater X = (F, W)^T
X = X + dt/6*(K1+2*K2+2*K3+K4);

% % theta differens og power supply til plot

Tdiff = X(idx(1),2)-X(idx(2), 2); % update the theta difference on current loop
Theta{end+1} = X(:,2); % add to theta difference array
P_supplied = A(2)*sin(Tdiff); % update the power supply on the current loop
Power(end+1) = P_supplied; % add to power supply array
Time(end+1) = t; % update time array
 
if t > time/2
    R(end + 1) = 1/N*real((sum(exp(1i*X(:,2)))));
    TimeR(end+1) = t - time/2;
end

end

end

