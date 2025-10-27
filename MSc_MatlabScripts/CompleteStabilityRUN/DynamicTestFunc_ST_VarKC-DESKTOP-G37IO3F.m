function [DATAarray, CELLarray] = DynamicTestFunc_ST_VarKC(MainType, SubType, epsilon)
% Perform the dynamic stability analysis of Type-"MainType"-"SubType" grid
% with epsilon being the treshold variance such that e/G^2 < epsilon where
% G is the number of simulated grids and e is the cummulated variance among
% the G grids. KC is the coupling strength factor multiplied upon all
% connection lines capacity (set to KC = 3 or 2).
% NOTE: MainType = 1,2,3,4, SubType = 'C','U','R'
%
% OUTPUT:
% - DATAarray: DATAarray{1,j,g} - g is the current grid (g=1..G), j=1 is
% the MagnitudeList of grid g, j = 2 is the connection matrix of grid g, j
% =3 is the power vector of grid g, j = 4 is the frequency of the
% magnitudes of magnitudelist, i.e. frequency(i) = MagnitudeList(i)/N_g
% where N_g is the number of nodes of grid g
%
% - CELLarray: CELLarray{1,i} - i = 1 is a magnitude vector
% (1+advance:advance:magMax), i=2 is an occurence vector MagFreq which contains all measured
% magnitudes among all grids of the grid-type and includes potential dublicates, that
% is,MagFreq(1:length(DATAarray{1,1,1}) for instance includes measured
% magnitudes of grid g=1 for all nodes (in order i=1...N_g)
% all simulated grids, i = 3 is a meanMagVec consisting of the magnitude
% means of each simulated grid (used to make normal distribution plot)

G = 1; % grids simulated
e = 0; % cummulated standard deviaion 
eps = epsilon;
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


while (e/G^2 >= eps || firstrun) && G < 5
    
    firstrun = 0;
    
    while ~bool
        [A,Pvec] = createTopology(MainType,SubType);
        N = length(Pvec);
        Omega0 = zeros(N,1);
        Theta0 = ones(N,1);


        % find and initiate steady state
        KR = 5*max(abs(Pvec)); %maxium trial capacity
        KL = 1; % minimum trial capacity
        KC = (KR+KL)/2; % initial coupling strength
        for j=1:7
            [bool, X] = simTraj(2,Time,dt,alpha,KC*A, Pvec', Omega0, Theta0);
            Omega0 = X(:,1);
            Theta0 = X(:,2);
            if bool
                KR = KC;
            else
                KL = KC;
            end
            KC = (KR+KL)/2;
        end
        
        KC = KR; % set KC to upper bound of interval to ensure synchrony
    end
    
    % MagnitudeList(i) contains the magnitude of node i in the current grid
    MagnitudeList = zeros(1,N);
    disp(G);
    
    for i=1:N
        
        P = Pvec;
        magR = round(max(abs(Pvec))/min(abs(Pvec)))+2;
        magL = 0;
        magnitude = magR;
        
        while abs(magL-magR)>1/3
            P(i) = magnitude*Pvec(i);
            bool = simTraj(2,Time,dt,alpha,KC*A,P',Omega0,Theta0);
            if bool
                magL = magnitude;
                magnitude = (magL+magR)/2;
            else
                magR = magnitude;
                magnitude = (magL+magR)/2;
            end
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
    DATAarray{1,2,G-1} = KC*A;
    DATAarray{1,3,G-1} = Pvec;
    DATAarray{1,4,G-1} = freq;
    
    NodesMeasured = NodesMeasured + N;
    
    % reset
    bool = 0;
        
end


% results for BOXPLOT and NORMAL DISTRIBUTION
% Mag contains vector of magnitudes (from 1+advance up till highest magnitude)
% MagFreq is a magnitude vector containing ALL measured magnitudes
% (including dublicates) 
% meanMagVec(i) contains the average number of magnitudes of grid i and we
% will use mean(meanMagVec) and std(meanMag)) to make a normal distribution
% NOTICE mean(meanMagVec) is normal distributed according to central limit theorem
Mag = 1+advance:advance:magMax;
MagFreq = zeros(1,NodesMeasured);
meanMagVec = zeros(1,G-1);
for i=1:G-1
    meanMagVec(i) = mean(DATAarray{1,1,i});
    idxMagFreq = find(MagFreq==0,1,'first');
    MagFreq(idxMagFreq:idxMagFreq + length(DATAarray{1,1,i})-1) = DATAarray{1,1,i};
end
std_meanMag = std(meanMagVec);
meanMag = mean(meanMagVec);

% create normaldistribution plot
figure()
x = 0:0.01:magMax;
norm = normpdf(x,meanMag, std_meanMag);
plot(x,norm)

% fill CELLarray for output
CELLarray = cell(1,3);
CELLarray{1,1} =Mag;
CELLarray{1,2} = MagFreq;
CELLarray{1,3} = meanMagVec;

end

