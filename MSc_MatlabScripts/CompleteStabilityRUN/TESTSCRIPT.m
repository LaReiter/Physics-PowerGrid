% testscript


G = 1; % grids simulated
e = 0; % cummulated standard deviaion 
eps = 1e-2;
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


while e/G^2 >= eps || firstrun
    
    firstrun = 0;
    
    while ~bool
        [A,Pvec] = createTopology_v2(1,'C',3);
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
        
        tic
        while bool
            magnitude = magnitude + advance;
            P(i) = magnitude*Pvec(i);
            bool = simTraj(2,Time,dt,alpha,A,P',Omega0,Theta0);
        end
        toc
        
%         P = Pvec;
%         magR = 60;
%         magL = 0;
%         magnitude = magR;
%         
%         tic
%         while abs(magL-magR)>1
%             P(i) = magnitude*Pvec(i);
%             bool = simTraj(2,Time,dt,alpha,A,P',Omega0,Theta0);
%             if bool
%                 magL = magnitude;
%                 magnitude = (magL+magR)/2;
%             else
%                 magR = magnitude;
%                 magnitude = (magL+magR)/2;
%             end
%         end
%         toc
        
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