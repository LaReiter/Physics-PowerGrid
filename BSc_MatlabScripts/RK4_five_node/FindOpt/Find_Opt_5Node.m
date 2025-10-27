clear all; close all;
%%% Simulation of the power grid kuramato model %%%
%%% Find optimal topology in the sense that it shows Braess Paradox

% ground parameters
dt = 1e-1; % timestep
totalTime = 400; % time duration
N0 = 5; % number of oscillators
found = 0;
eps = 0.2;
pmaxMax = 5;

while found == 0;
    
    Rinf = 0;
    while Rinf < eps
        pmax = randi(pmaxMax);
        
        % Create random N0 node coupling matrix
        A = CreateAdj(N0); 

        % initial conditions
        W0 = zeros(N0,1);
        T0 = zeros(N0,1);
        P0 = ceil(rand(1,N0-round(N0/2))*randi(2));
        s1 = sum(P0);
        PG = repmat(s1/N0, 1, round(N0/2));
        P0 = transpose([-P0 PG]);
        alpha = rand(1); % dissipation parameter

        [Kc, Rinf]  = findKc(totalTime, dt, alpha, A, P0, W0, T0);
        
    end
    
    disp('check1')

    idx1 = find(triu(A));
    % check Braess paradox on links
    Links = []; % BP visible
    for k=idx1.'
        A0 = A;
        A0(k) = 2*A0(k);
        A0 = triu(A0) + triu(A0).';
        [Kc_Links, Rinf_Links] = findKc(totalTime, dt, alpha, A0, P0, W0, T0);
        if (Rinf_Links > eps) && (Kc_Links > Kc);
            Links(end+1) = k
        end
    end

    disp('check2')
    if sum(Links) == 0
        continue
    end;

    idx2 = find(~(triu(A)+tril(ones(N0,N0))));
    % check Braess Paradox on new links
    PotLinks = []; % BP visible
    for j=idx2.'
        A0 = A;
        A0(j) = pmax;
        A0 = triu(A0) + triu(A0).';
        [Kc_PotLinks, Rinf_PotLinks] = findKc(totalTime, dt, alpha, A0, P0, W0, T0);
        if (Rinf_PotLinks > eps) && (Kc_PotLinks > Kc)
            PotLinks(end+1) = j
        end
    end
    
    
    if sum(PotLinks) == 0
        continue
    else
        found = 1; % Else topology found    
    end
    
end
