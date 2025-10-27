clc; clear all; close all;
% Lars Reiter Nielsen
% Kuramoto model: en forbruger og en generator
% Sammenligning af stabilt fixpunkt og instabilt fixpunkt
% instabile linjer konvergerer mod stabilt fixpunkt

% Grundparametre
dt = 1e-1; % BØR VÆLGES SMART
alpha = 0.1; % dissipations parameter
N = 2; % antal oscillatorer
time = 400;
P = [1; -1];
dP = 2;
idx = [1,2];
pmax = 3;   % stabiliserer dp = -3.5 ved pmax = 3.138
A = pmax*[[0 1]; [1 0]];




% Definer tilstands vektor X = (W,Theta)^T således at systemet hedder 
% dX/dt = F(X)
unitV = ones(N,1);

% % definer F(X)
f = @(V, P) [(P - alpha*V(:,1) + (A.*sin(transpose(repmat(V(:,2),1,N)) - repmat(V(:,2),1,N)))*unitV) V(:,1)];

cellTheta = {};
hold on
for i=1:2
    % tids og faseforskels vektor til plotning
    DTheta = [];
    Time = [];

    % frekvenser: theta_i er fasen af oscillator i og w_i er fasehastigheden
    W = [0;0];
    if i == 1
        T = [1/2*asin(dP/(2*pmax));-1/2*asin(dP/(2*pmax))];
    else
        T = [1/2*(pi-asin(dP/(2*pmax)));-1/2*(pi-asin(dP/(2*pmax)))];
    end
    
    X = [W T];

    for t = 0:dt:time

    Ps = P;
    
    % Pertubation
    if (t < 42) && (40 < t)
        Ps = P + [0;-1e-1];
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

    % time
    Time(end+1) = t; % opdater tiden

    end
    
    DTheta2 = [];
    Time = [];

    % frekvenser: theta_i er fasen af oscillator i og w_i er fasehastigheden
    W = [0;0];
    if i == 1
        T = [1/2*asin(dP/(2*pmax));-1/2*asin(dP/(2*pmax))];
    else
        T = [1/2*(pi-asin(dP/(2*pmax)));-1/2*(pi-asin(dP/(2*pmax)))]
    end
    X = [W T];

    for t = 0:dt:time

    Ps = P;
    
    % Pertubation
    if (t < 202) && (t > 200);
        Ps = P + [0;-2];
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
    DTheta2(end+1) = Tdiff; % tilføj til listen

    % time
    Time(end+1) = t; % opdater tiden

    end
    
    cellTheta{end+1} = [DTheta; DTheta2];
end

[hAx, hLine1, hLine2] = plotyy(Time, cell2mat(cellTheta(1)), Time, cell2mat(cellTheta(2)));
title('Stability against pertubation for fixpoint F_1 and F_2');
xlabel('time (a.u.)');
ylabel(hAx(1), '\Delta \theta (F1 fixpoint)')
hold off

