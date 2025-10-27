clc; clear all; close all;
% Lars Reiter Nielsen
% Kuramoto model: en forbruger og en generator
% Genkendelse af fixpunkter og stabilitetsanalyse

% Grundparametre
dt = 1e-1; % BØR VÆLGES SMART
alpha = 0.1; % dissipations parameter
N = 2; % antal oscillatorer
time = 200;
P = [1; -1];
dP = [0; -3.4];
idx = [1,2];
pmax = 3;   % stabiliserer dp = -3.5 ved pmax = 3.138
A = pmax*[[0 1]; [1 0]];

% frekvenser: theta_i er fasen af oscillator i og w_i er fasehastigheden
T = [29.3596;29.0202];
W = [-0.001;-0.001];

% Definer tilstands vektor X = (W,Theta)^T således at systemet hedder 
% dX/dt = F(X)
X = [W T];
unitV = ones(N,1);

% % definer F(X)
f = @(V, P) [(P - alpha*V(:,1) + (A.*sin(transpose(repmat(V(:,2),1,N)) - repmat(V(:,2),1,N)))*unitV) V(:,1)];


% tids og faseforskels vektor til plotning
DTheta = [];
DOmega = [];
Time = [];
Power = [];

for t = 0:dt:time

Ps = P;

% Pertubation
if (t < 42) && (40 < t)
    Ps = P + dP;
end
    
% opdater approksimationer til f(X)
a = f(X, Ps);
b = f(X+dt/2*a, Ps);
c = f(X+dt/2*b, Ps);
d = f(X+dt*c, Ps);

% opdater X = (F, W)^T
X = X + dt/6*(a+2*b+2*c+d);

% theta differens
Tdiff = X(idx(1),2)-X(idx(2), 2); % opdater til relevante theta differens 
DTheta(end+1) = Tdiff; % tilføj til listen

% power supply
P_supplied = A(2)*sin(Tdiff); % opdater til relevante power supply
Power(end+1) = P_supplied; % tilføj til listen

% omega differens
Odiff = X(idx(1),1)-X(idx(2),1); % opdater til relevant omega differens
DOmega(end+1) = Odiff;

% time
Time(end+1) = t; % opdater tiden

end

hold on
[hAx, hLine1, hLine2] = plotyy(Time, DTheta, Time, DOmega);
title('No resynchronization. Maximimum capacity exceeded.');
xlabel('time (a.u.)');
ylabel(hAx(1), '\theta_1 - \theta_2');
ylabel(hAx(2), '\omega_1 - \omega_2');
hold off