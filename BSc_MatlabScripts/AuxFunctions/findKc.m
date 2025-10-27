function [ Kc, Rinf ] = findKc( time, dt, alpha, A, P0, W0, T0)
% Finds and returns the critical coupling strength KC based on the topology
% and simulation setup.
% time: Total simulation time
% dt: timestep
% alpha: dissipation parameter
% mat: coupling (or adjacency) matrix
% p0: initial power vector
% w0: initial phase velocity vector
% t0: initial phase vector

eps = 1e-3;
pmax = 20;
pmin = 0;
pdiff = pmax;
Kc = 0;

while (pdiff > eps)
    
    pp = Kc;  % save previous coupling point
    Kc = (pmax+pmin)/2; % new coupling point
    pdiff = abs(Kc - pp); % difference from previous and current coupl str

    [Rinf] = NnodeSim_CritVal(time,dt, alpha, Kc*A,P0, W0, T0);
    
    if (Rinf < 0.05) && (abs(Kc) < eps)
        Kc = -1;  % intercept and conclude no steady state
        return;
    end
        
    if Rinf > 0.05 % steady state reached
        pmax = Kc;
    else
        pmin = Kc;
    end

end



end

