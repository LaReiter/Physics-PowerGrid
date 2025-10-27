clc; clear all; clf;
%%% Simulating dynamic stability

G = 1; % grids simulated
e = 0; % cummulated standard deviaion 
epsilon = 1e-2;
advance = 1/4;
Time = 100;
dt=1e-1;
alpha=0.1;
magMax = 0;

NodesMeasured = 0;
bool = 0;

firstrun = 1;

% DATAarray{i,j,k} with i=1, j = 1,2,3,4 and k=1,2,3,...,M and j=1 is
% MagnitudeList of grid simulation k, j=2 is A of grid simulation k, and
% j=3 is Pvec of simulation k, j = 4 is the frequency of the MagnitudeList
M = 200;
DATAarray = cell(1,3,M);


while e/G^2 >= epsilon || firstrun
    
    firstrun = 0;
    
    while ~bool
        [A,Pvec] = createTopology_v2(4,'C',3);
        N = length(Pvec);
        Omega0 = zeros(N,1);
        Theta0 = ones(N,1);


        % find and initiate steady state
        [bool, X] = simTraj(2,Time,dt,alpha,A, Pvec', Omega0, Theta0);
        Omega0 = X(:,1);
        Theta0 = X(:,2);
    end
    
    % MagnitudeList(i) contains the magnitude of node i in the current grid
    MagnitudeList = zeros(1,N);
    disp(G);
    
    for i=1:N
        
        P = Pvec;
        magnitude = 1;
        
        while bool
            magnitude = magnitude + advance;
            P(i) = magnitude*Pvec(i);
            bool = simTraj(2,Time,dt,alpha,A,P',Omega0,Theta0);
        end
        
        magMax = max(magMax,magnitude); % update maximum magnitude
        
        MagnitudeList(i) = magnitude;
        % reset
        bool = 1;
        
                    
    end
    
    % update standard error and advance grid number
    e = e + var(MagnitudeList);
    G = G + 1;
    
    % freq(i) = "occurences of MagnitudeList(i)"/N
    freq = zeros(size(MagnitudeList));
    for i = 1:length(MagnitudeList)
        freq(i) = sum(MagnitudeList==MagnitudeList(i));
    end
    freq = 1/N*freq;
    
    disp(e/G^2);
    % save data
    DATAarray{1,1,G-1} = MagnitudeList;
    DATAarray{1,2,G-1} = A;
    DATAarray{1,3,G-1} = Pvec;
    DATAarray{1,4,G-1} = freq;
    
    NodesMeasured = NodesMeasured + N;
    
    % reset
    bool = 0;
        
end


% results for BOXPLOT and NORMAL DISTRIBUTION
% Mag contains vector of magnitudes (from 1+advance up till highest magnitude)
% MagFreq is a magnitude vector with MagFreq(i) = 1 + i*advance (magnitude
% strength) uptill magnitude maximum divideded by advance
% meanMagVec(i) contains the average number of magnitudes of grid i and we
% will use mean(meanMagVec) and std(meanMag)) to make a normal distribution
% NOTICE mean(meanMagVec) is normal distributed according to central limit theorem
Mag = 1+advance:advance:magMax;
MagFreq = zeros(1,magMax/advance);
meanMagVec = zeros(1,G-1);
for i=1:G-1
    [C,Cindices] = unique(DATAarray{1,1,i});
    meanMagVec(i) = mean(DATAarray{1,1,i});
    for j=1:length(C)
        MagFreq((C(j)-1)/advance) = MagFreq((C(j)-1)/advance) + sum(DATAarray{1,1,i}==C(j)); % update magnitude vector (add number of occurences of magnitude C(j) for grid i)
    end
end
std_meanMag = std(meanMagVec);
meanMag = mean(meanMagVec);

% create normaldistribution plot
figure()
x = 0:0.01:magMax;
norm = normpdf(x,meanMag, std_meanMag);
plot(x,norm)

% freq = 1/NodesMeasured*MagnitudeList;
% mag = 1+advance:advance:maxPower;
% 
% bar(mag,freq)
% xlim([1+advance maxPower])