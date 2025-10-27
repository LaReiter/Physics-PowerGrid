function [A,Pvec] = createTopology_v2( Type, Subtype, CS )
% createTopology: Creates a topology of Type-'Type' and Subtype 'Subtype'
% EDITED v2: Added variable connection strength to adjacency matrix. Select
% connected components and measure max-power Pmax (numerically), let K =
% CS*Pmax (price-optimal, real-world inspired)
% Lars Reiter Nielsen

% average degree (of standard network)
d = 2.7;
th_min = 0.4; % lowerbound threshold for huge oscillator distribution
th_max = 0.8; % upperbound threshold for huge oscillator distribution

% Define total power and power vectors
P = 25;
Pc = [1/5,1/2, 1, 2, 6];
Pg = [1/5,33/20, 31/10, 91/20, 6];

% Define GBP and CBP
delta = zeros(1,5);
mu = zeros(1,5);

% normal delta and mu vector and tolerance vector
deltamean = [0.5,0.2,0.15,0.1,0.05];
deltadev = [0.2,0.1,0.1,0.02];
mumean = [0.5,0.2,0.15,0.1,0.05];
mudev = [0.2,0.1,0.1,0.02];

% Choose CBP and GBP according to type
switch Type
    case 1
        for i=1:4
            if i > 1
                corr_delta = delta(i-1)-(deltamean(i-1)-corr_delta);
                corr_mu = mu(i-1) - (mumean(i-1)-corr_mu);
            else
                corr_delta = 0;
                corr_mu = 0;
            end
            lower_delta = deltamean(i)-corr_delta - 1/2*deltadev(i);
            higher_delta = deltamean(i)-corr_delta + 1/2*deltadev(i);
            delta(i) = (higher_delta-lower_delta).*rand(1) + lower_delta;
            
            lower_mu = mumean(i)-corr_mu - 1/2*mudev(i);
            higher_mu = mumean(i)-corr_mu + 1/2*mudev(i);
            mu(i) = (higher_mu-lower_mu).*rand(1) + lower_mu;
            
        end

        delta(5) = 1-sum(delta);
        mu(5) = 1-sum(mu);        
    case 2
        for i=1:4
            if i > 1
                corr_delta = delta(7-i)-(deltamean(i-1)-corr_delta);
                corr_mu = mu(i-1) - (mumean(i-1)-corr_mu);
            else
                corr_delta = 0;
                corr_mu = 0;
            end
            lower_delta = deltamean(i)-corr_delta - 1/2*deltadev(i);
            higher_delta = deltamean(i)-corr_delta + 1/2*deltadev(i);
            delta(6-i) = (higher_delta-lower_delta).*rand(1) + lower_delta;
            
            lower_mu = mumean(i)-corr_mu - 1/2*mudev(i);
            higher_mu = mumean(i)-corr_mu + 1/2*mudev(i);
            mu(i) = (higher_mu-lower_mu).*rand(1) + lower_mu;
            
        end

        delta(1) = 1-sum(delta);
        mu(5) = 1-sum(mu);
    case 3
        for i=1:4
            if i > 1
                corr_delta = delta(i-1)-(deltamean(i-1)-corr_delta);
                corr_mu = mu(7-i) - (mumean(i-1)-corr_mu);
            else
                corr_delta = 0;
                corr_mu = 0;
            end
            lower_delta = deltamean(i)-corr_delta - 1/2*deltadev(i);
            higher_delta = deltamean(i)-corr_delta + 1/2*deltadev(i);
            delta(i) = (higher_delta-lower_delta).*rand(1) + lower_delta;
            
            lower_mu = mumean(i)-corr_mu - 1/2*mudev(i);
            higher_mu = mumean(i)-corr_mu + 1/2*mudev(i);
            mu(6-i) = (higher_mu-lower_mu).*rand(1) + lower_mu;
            
        end

        delta(5) = 1-sum(delta);
        mu(1) = 1-sum(mu);        
    case 4
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
end


% Define generator and consumer size vector (entry = small to large)
Nc = floor((mu ./ Pc)*P); 
Ng = floor((delta ./ Pg)*P);

% Adjust with rescale parameter
beta = Pc*transpose(Nc)/(Pg*transpose(Ng));
Pg = beta*Pg;

% total edges to be distributed
N = sum(Nc) + sum(Ng);
E = floor(N*d/2);

% define power vector
% Notice Pvec is sorted, such first entries contains tiny consumers, then
% tiny generators, then small consumers, then small generators and so
% forward. Medium generators can for instance be found as
% Pvec(Nc(1)+Ng(1)+Nc(2)+Nc(2)+Nc(3):Nc(1)+Ng(1)+Nc(2)+Nc(2)+Nc(3)+Ng(3))
Pvec = zeros(1,N);
idx = 1;
for i=1:5
    Pvec(idx:idx+Nc(i)-1) = -Pc(i);
    Pvec(idx+Nc(i):idx + Nc(i) + Ng(i) -1) = Pg(i);
    idx = idx + Nc(i) + Ng(i);

end


% construct frame for adjacency matrix
% once again the first entries contain tiny consumers, then tiny generators
% then small consumers and so forward (see power vec).
A = createMinSpanTree(N,N/2);

E0 = E - (N-1);

 switch Subtype
    % TODO: How to handle case Nl == 0? Is it already handled?
    % Check "Speciale" section "Characterising quasi-isomorphic topologies"
    % for more information
    case 'C'
        
        A = zeros(N,N);
        while ~isConnected(A)
            A = createBiConGraph(N,E); 
        end


    case 'U'
        
        edgelist = sum(A);        
        curNode = 1;
        
        while E0>0
            [minel,minidx] = min([edgelist(1:curNode-1),edgelist(curNode+1:end)]);
            
            % adjust indices
            if minidx >= curNode
                minidx = minidx + 1;
            end
            
            A(curNode,minidx) = 1;
            A(minidx,curNode)=1;
            
            edgelist(curNode)=edgelist(curNode)+1;
            edgelist(minidx)=edgelist(minidx)+1;
            
            E0 = E0-1;
            
            [minel,curNode] = min(edgelist);
        end

        
    case 'R'
        
        sample2 = find(triu(ones(N,N),1)>0);

        
        % remove diagonal elements from sample
        Erandom = datasample(sample2, E0, 'Replace',false);
        A(Erandom) = 1;
                
                

end

A = A + A';
% A(A>1) = 1;

% Implement variable connection strength
for i = 1:N
    indices = find(A(i,:)>0);
    for j=1:numel(indices)
        A(i,indices(j)) = CS*max(abs(Pvec(i)),abs(Pvec(indices(j))));
    end
end

end

