% testscript
% Define GBP and CBP
delta = zeros(1,5);
mu = zeros(1,5);

% normal delta and mu vector and tolerance vector
deltamean = [0.5,0.2,0.15,0.1,0.05];
deltadev = [0.2,0.1,0.1,0.02];
mumean = [0.5,0.2,0.15,0.1,0.05];
mudev = [0.2,0.1,0.1,0.02];


        for i=1:4
            if i > 1
                corr_delta = delta(7-i)-(deltamean(i-1)-corr_delta);
                corr_mu = mu(7-i) - (mumean(i-1)-corr_mu);
            else
                corr_delta = 0;
                corr_mu = 0;
            end
            lower_delta = deltamean(i)-corr_delta - 1/2*deltadev(i);
            higher_delta = deltamean(i)-corr_delta + 1/2*deltadev(i);
            delta(6-i) = (higher_delta-lower_delta).*rand(1) + lower_delta;
            
            lower_mu = mumean(i)-corr_mu - 1/2*mudev(i);
            higher_mu = mumean(i)-corr_mu + 1/2*mudev(i);
            mu(6-i) = (higher_mu-lower_mu).*rand(1) + lower_mu;
            
        end

        delta(1) = 1-sum(delta);
        mu(1) = 1-sum(mu);