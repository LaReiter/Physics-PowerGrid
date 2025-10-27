clc; clear all; clf;
%%% Simulating dynamic stability

grids = 50;
advance = 1/4;
Time = 100;
dt=1e-1;
alpha=0.1;
maxPower = 10;

MagnitudeList = zeros(1,(maxPower-1)/advance);
NodesMeasured = 0;
bool = 0;


for grid=1:grids
    
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
    
    disp(grid);
    
    for i=1:N
        
        
        P = Pvec;
        magnitude = 1;
        
        while bool
            magnitude = magnitude + advance;
            P(i) = magnitude*Pvec(i);
            bool = simTraj(2,Time,dt,alpha,A,P',Omega0,Theta0);
        end
        
        idx = round((magnitude-1-advance)/advance);
        MagnitudeList(1:idx) = MagnitudeList(1:idx) + 1;
        bool = 1;
        
                    
    end
    
    NodesMeasured = NodesMeasured + N;
        
end


freq = 1/NodesMeasured*MagnitudeList;
mag = 1+advance:advance:maxPower;

bar(mag,freq)
xlim([1+advance maxPower])