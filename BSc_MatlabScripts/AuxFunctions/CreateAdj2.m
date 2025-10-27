    function [ A ] = CreateAdj2( N,pmax, V, M)
% Code by Lars Reiter Nielsen @ Student at Copenhagen University
% N: Number of nodes
% V: Number of edges
% M: Non-complete base adjacency matrix for fixed edges
% Creates a random N node adjacency matrix (undirected - connected - graph) 
% with V edges. If M matrix is given as input, entries in M will apply to A


switch nargin
    case 2 % random N node topology
        A = zeros(N,N); % define base matrix
        for i=1:N-1
            A(i,i+randi(N-i)) = pmax; % add N-1 edges connecting every node
        end
        upA = triu(round(rand(N))*pmax,1);
        A = A + A.' + upA + upA.';
        A(A>pmax) = pmax; % adjust for double edges
        
    case 3 % N node topology with V edges
        A = zeros(N,N); % define base matrix
        for i=1:N-1
            A(i,i+randi(N-i)) = pmax; % add N-1 edges connecting every node
        end
        curV = length(find(A)); % number of current edges
        missV = V-curV; % determine missing edges
        potV = find(~(A+tril(ones(N,N)))); % find potential edges
        if missV > 0 % missing edges, add
            newV = datasample(potV, missV);
            A(newV) = pmax;
            A = A + A.';
        else % missV = 0
            A = A + A.';
        end
        
    case 4 % N node topology with V edges based on M with pmax crit val
        M0 = [M zeros(30, N-30); zeros(N-30,N)];
        A = triu(M0,1); % define base matrix
        for i=1:N-1
            if isempty(find(A(i,:), 1))
                A(i,i+randi(N-i)) = pmax; % add N-1 edges connecting every node
            end
        end
        curV = length(find(A)); % number of current edges
        missV = V-curV; % determine missing edges
        potV = find(~(A+tril(ones(N,N)))); % find potential edges
        if missV > 0 % missing edges, add
            newV = datasample(potV, missV);
            A(newV) = pmax;
            A = A + A.';
        else % missV = 0
            A = A + A.';
        end

end

