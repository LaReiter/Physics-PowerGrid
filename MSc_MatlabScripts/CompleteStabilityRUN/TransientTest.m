clc; clear all; clf;
%%% Simulating transient (basin) stability

bool = 0;

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