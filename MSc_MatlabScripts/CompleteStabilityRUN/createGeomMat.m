function [ A ] = createGeomMat( N, shape)
% 
% N: Nodes of network
% Shape: Determines the shape of the network. Choose between:
% 'Cyclic': Network where every node has at least degree 2. No inner trees.
% 'Dead': Together with 'noLeaves', creates a network where 'noLeaves'
% defines the number of (random) nodes which are leaves of the network and
% thus are attached to exactly one node.
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

