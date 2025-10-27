function [A,Pvec] = createTopology( Type, Subtype )
% createTopology: Creates a topology of Type-'Type' and Subtype 'Subtype'
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




% OLD code for another implementation of the CBP and GBP
% switch Type
%     case 1
%         delta(1) = th_min+(th_max-th_min)*rand();
%         mu(1) = th_min+(th_max-th_min)*rand();
%         for i = 2:5
%             delta(i) = (1-sum(delta))*rand();
%             mu(i) = (1-sum(mu))*rand();
%         end
%     case 2
%         delta(5) = th_min+(th_max-th_min)*rand();
%         mu(1) = th_min+(th_max-th_min)*rand();
%         for i = 2:5
%             delta(6-i) = (1-sum(delta))*rand();
%             mu(i) = (1-sum(mu))*rand();
%         end
%     case 3
%         delta(1) = th_min+(th_max-th_min)*rand();
%         mu(5) = th_min+(th_max-th_min)*rand();
%         for i = 2:5
%             delta(i) = (1-sum(delta))*rand();
%             mu(6-i) = (1-sum(mu))*rand();
%         end
%     case 4
%         delta(5) = th_min+(th_max-th_min)*rand();
%         mu(5) = th_min+(th_max-th_min)*rand();
%         for i = 2:5
%             delta(6-i) = (1-sum(delta))*rand();
%             mu(6-i) = (1-sum(mu))*rand();
%         end
% end

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
%     idx = find(Pvec==0,1);
%     if Nc(i) ~= 0
%         Pvec(idx:idx+Nc(i)-1) = Pc(i);
%         if Ng(i) ~= 0
%             Pvec(idx+Nc(i):idx+Nc(i)+Ng(i)) = Pg(i);
%         end
%     elseif Ng(i)~= 0
%         Pvec(idx:idx+Ng(i)-1) = Pg(i);
%     end
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
    case "C"
        
        A = zeros(N,N);
        while ~isConnected(A)
            A = createBiConGraph(N,E); 
        end
        
%         % size factor for larger oscillators
%         sfac = 3;
%         
%         Nl = Nc(4) + Nc(5) + Ng(4) + Ng(5);
%         Ns = N - Nl;       
%         
%         % scheme determines if we switch to small - to - large connection scheme
%         scheme = [0,0];
%         
%         % smaller oscillators indices (tiny, small,medium)
%         Sindices = zeros(1,Ns^2);
%         for n = 0:(Ns-1)
%             Sindices(1+Ns*n:Ns+Ns*n) = (1+n*N):(N*n+Ns);
%         end
%         
%         % larger oscillators indices (large, huge)
%         Lindices = zeros(1,Nl^2);
%         for n = 0:(Nl - 1)
%             Lindices(1 + Nl*n:Nl+Nl*n) = (N+1)*Ns+1+n*N:(Ns+(n+1))*N;
%         end
%         
%         Mixindices = setdiff(setdiff(1:N^2,Lindices),Sindices);
%         
%         % remove diagonal elements
%         diagindices = 1:(N+1):N^2;
%         Sindices = setdiff(Sindices,diagindices);
%         Lindices = setdiff(Lindices, diagindices);
%         Mixindices = setdiff(Mixindices,diagindices);
%         
%         % find indices of entries filled
%         FilledEntries = find(A>0);
%         
%         % Adjust and remove already filled entry indices
%         Sindices = setdiff(Sindices, FilledEntries);
%         Lindices = setdiff(Lindices, FilledEntries);
%         
%         count = floor(E0/2);
%        
%         if Nl ~= 0
%             while count > 0
%                El = sum(sum(A(:,N-Nl+1:N)));
% %                El = sum(sum(A(:,end:-1:N-Nl+1)));
%                Es = sum(sum(A(:,1:N-Nl)));
%                
%                if El/Nl > sfac*Es/Ns && scheme(2) < 1
%                    
%                    possibleNodes = setdiff(Sindices,FilledEntries);
%                    
%                    if isempty(possibleNodes)
%                        scheme(2) = 1;
%                    else                  
%                        % add connection between two small
%                        node = datasample(possibleNodes,1);
%                        A(node) = 1;
%                        A = A + A';
%                        A(A>1) = 1;
%                    end
%                    
%                    FilledEntries = find(A>0);
%                    scheme(1) = 1;
%                elseif El/Nl <= sfac*Es/Ns && scheme(1) < 1
%                    
%                    possibleNodes = setdiff(Lindices,FilledEntries);
%                    
%                    if isempty(possibleNodes)
%                        scheme(1) = 1;
%                    else
%                        % add connection between two large
%                        node = datasample(possibleNodes,1);
%                        A(node) = 1;
%                        A = A + A';
%                        A(A>1) = 1;
%                    end
%                    
%                    FilledEntries = find(A>0);
%                    
%                    scheme(2) = 1;
%                else
%                    % add connection between small and large
%                    node = datasample(setdiff(Mixindices,FilledEntries),1);
%                    A(node) = 1;
%                    A = A + A';
%                    A(A>1) = 1;
%                    
%                    FilledEntries = find(A>0);
%                end
%                       
% 
%                count = count - 1;
%             end
%         else
%         end
        

    case "U"
        
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
        
        
%         auxMat = createMinSpanTree(N,N);
% 
%         size = E; 
%         i = 1;
%         
%         while size > 0
% %            filledEntries = find(auxMat(i,:)>0);
% %            sampleSize = max(1,min(round(E/N)-length(filledEntries),size));
% % 
% %            auxMat(i,datasample(setdiff(1:N,filledEntries),sampleSize, 'Replace',false)) = 1; 
% %            auxMat = auxMat + auxMat';
% 
%            freeEntries = find(auxMat(i,:)==0);
%            auxMat(i,datasample(freeEntries,1, 'Replace',false)) = 1; 
%            auxMat = auxMat + auxMat';
% 
%            
%            size = size - 1;
%            i = max(mod(i+1,N+1),1);
%         end
%         
%         A = auxMat;
        
    case "R"
        
        sample2 = find(triu(ones(N,N),1)>0);
        
%         % Diagonal element list
%         DiagElem2 = 1:N;
%         DiagElem2 = DiagElem2*(N+1)-N;
        
        % remove diagonal elements from sample
        Erandom = datasample(sample2, E0, 'Replace',false);
        A(Erandom) = 1;
                
                

end

A = A + A';
A(A>1) = 1;


% % %% show average degree
% % d = sum(sum(A))/N;
% % disp(d)

% % show El/Nl and Es/Ns
% disp(El/Nl)
% disp(Es/Ns)
end

