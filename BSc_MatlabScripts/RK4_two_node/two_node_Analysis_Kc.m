clc; clear all; close all;
% Lars Reiter Nielsen
% Kuramoto model: en forbruger og en generator
% Bestem optimal koblingsstyrke for de to stabile fixpunkter
% hhv. T_1 = (0, arcsin(\Delta P)) og T_2 = (0, pi - arcsin(\Delta P))

% Grundparametre
dt = 1e-1; % BØR VÆLGES SMART
N = 2; % antal oscillatorer
time = 1500;
A = [[0 1]; [1 0]];


% Definer tilstands vektor X = (W,Theta)^T således at systemet hedder 
% dX/dt = F(X)
unitV = ones(N,1);

% % definer F(X)
f = @(V, P,alpha,k) [(P - alpha*V(:,1) + (k*A.*sin(transpose(repmat(V(:,2),1,N)) - repmat(V(:,2),1,N)))*unitV) V(:,1)];

hold on
% power difference
DP = [];

linewidth = 0;
% iterate over consumption
for pc = 1:0.5:2
    
    P = [-pc; pc];

    % forholds og gennemsnitlig fase ordens parameter liste
    Diff = [];
    Rinf = [];
    Time = [];
    alpha = 0.1;
    
    % max kobling og forventet kritisk punkt
    pmax = 0;
    expFixP = (P(2)-P(1));
    
    while pmax < 5

        % forøg pmax
        if (pmax > expFixP-0.3) && (pmax < expFixP+0.3)
            pmax = pmax + 0.05;
        else
            pmax = pmax + 0.2;
        end

        % fixpunkt
        T = [asin((P(2)-P(1))/pmax); asin((P(2)-P(1))/pmax)];
        W = [0;0];

        X = [W T];
        Diff(end+1) = pmax;
        R = []; % fase ordens parameter liste

        for t = 0:dt:time  

            Ps = P;

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

    
    plot(Diff, Rinf,'LineWidth',2);
    DP(end+1) = expFixP;
    
    
end

title('Critical coupling strength')
xlabel('2K')
ylabel('r_{\infty}            ','rot',0)
axis([1 5 -0.01 1])
legend(['\Delta P =' num2str(DP(1))],['\Delta P =' num2str(DP(2))],['\Delta P =' num2str(DP(3))])
hold off