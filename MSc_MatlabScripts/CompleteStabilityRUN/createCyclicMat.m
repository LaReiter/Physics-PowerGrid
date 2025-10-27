function [ A ] = createCyclicMat( N, E)
% 
% N: Nodes of network
% E: Number of edges to distribute
% create a cyclic network (with no dead ends) with N node and E edges
% Code by Lars Reiter Nielsen @ Student at Copenhagen University

A = zeros(N,N); % define base matrix

visited = linspace(1,N,N); % number of visited states must equal zero

if strcmp(shape, 'Cyclic')

    now = randi(N,1);
    prev = now;
    next = now;

    A(1,now) = 1;

    while sum(visited) > 0

        while next == prev
            next = randi(N,1);
        end

        A(now,next) = 1;

        visited(now)=0;

        prev = now;
        now = next;
        
        
    end

    while next == prev
        next = randi(N,1);
    end

    A(now, next)=1;

    disp(next)
    disp(now)
    disp(prev)
    
    A = A - diag(diag(A));  
    A = A+A';
    A(A>1) = 1;       


end


end

