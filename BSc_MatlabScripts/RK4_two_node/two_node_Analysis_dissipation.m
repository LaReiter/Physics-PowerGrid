clc; clear all; close all;
% Lars Reiter Nielsen
% Kuramoto model: en forbruger og en generator
% Script til bestemmelse af kritisk koblingsstyrke med pertubation

% Grundparametre
dt = 1e-1;
N = 2; % antal oscillatorer
time = 1500;
P = [1; -1];
idx = [1,2];
% pmax = 2; 
dp = -3.5;
A = [[0 1]; [1 0]];

% frekvenser: theta_i er fasen af oscillator i og w_i er fasehastigheden
T = [-29.3596;-29.0202];
W = [-0.0001;0.0001];
alpha = 0.1;
% T = [0; 0];
% W = [-0.0001;0.0001];

% Definer tilstands vektor X = (W,Theta)^T således at systemet hedder 
% dX/dt = F(X)
unitV = ones(N,1);

% % definer F(X)
f = @(V, P,alpha,k) [(P - alpha*V(:,1) + (k*A.*sin(transpose(repmat(V(:,2),1,N)) - repmat(V(:,2),1,N)))*unitV) V(:,1)];


% Kritisk koblingstyrke liste og gennemsnitlig fase ordens parameter liste
KcList = [];
Rinf = [];
Time = [];
for pmax = 1:0.2:5
X = [W T];
KcList(end+1) = pmax;
R = []; % fase ordens parameter liste

    for t = 0:dt:time  

    Ps = P;

    % Pertubation
    if (t < 42) && (40 < t)
        Ps = P + [dp;0];
    end
    
    % opdater approksimationer til f(X)
    a = f(X, Ps,alpha,pmax);
    b = f(X+dt/2*a, Ps,alpha,pmax);
    c = f(X+dt/2*b, Ps,alpha,pmax);
    d = f(X+dt*c, Ps,alpha,pmax);

    % opdater X = (F, W)^T
    X = X + dt/6*(a+2*b+2*c+d);

    if t > time/2
        R(end + 1) = 1/2*(exp(1i*X(1,2))+exp(1i*X(2,2)));
        Time(end+1) = t;
    end

    end
    disp(pmax)
Rinf(end+1) = real(sum(R))/length(R);
end

hold on
x_p = [3.5,3.5];
y_p = [0,1];
plot(x_p,y_p,'--')
plot(KcList, Rinf,'LineWidth',2);
title('Resynchronization through strong coupling')
xlabel('P_{max}')
ylabel('r_{\infty}            ','rot',0)
legend('Pertubation strength','Location', 'northwest')
axis([1 5 -0.01 1])
hold off