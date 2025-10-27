function [ Tarray ] = TransientTestFunc( MainType, SubType, epsilon, KC )
% Perform the transient stability analysis of Type-"MainType"-"SubType" grid
% with epsilon being the treshold variance.

    dist = 'Uniform'; % make distribution and generation equal for both consumers and generators
    Time = 100; % runtime
    dt = 1e-1; % timestep
    alpha = 0.1; % dissipation parameter
    Trials = 50; % number of pertubation trials in onenodeBasinStability
    G = 1; % initiate grid counter
    error = 100; % threshold error
    varSg = 0; % initiate sum of standard errors for grid g
    firstRun = 0;
    
    Tarray = cell(1,3,200); % cell{1,1,g} contains S or grid g, cell{1,2,g} contains A of grid g, cell{1,3,g} contains Pvec of grid g

    while ((1/G)*sqrt(error) > epsilon) && G <= 2
        
        disp(G)

        bool = 0;
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

        S = onenodeBasinStability(KC*A,Theta0,Omega0,Pvec',alpha,dt,Time,Trials,dist);
        Ng = length(S);
        varSg = sum(S.*(1-S))/Trials;

        if G==1
            error = varSg/(Ng^2);
        else
            error = error + varSg/(Ng^2);
        end
        

        if firstRun
            G = G+1; % advance grid counter
        end
        
        firstRun = 1; % firstrun elapsed
        
        % update cellarray
        Tarray{1,1,G} = S;
        Tarray{1,2,G} = KC*A;
        Tarray{1,3,G} = Pvec;
        
        disp((1/G)*sqrt(error))

    end

end

