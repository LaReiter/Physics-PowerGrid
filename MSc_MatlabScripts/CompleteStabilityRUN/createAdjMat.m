function [ A ] = createAdjMat( N, d )
%createAdjMat creates a N x N connection matrix with average degree d
%   Input parameters:
%   N:  Dimension of the square (symmetric) connection/adjacency matrix
%   d:  The average degree of the connection matrix (average number of
%   edges protruding from each node)
% Code by Lars Reiter Nielsen @ Student at Copenhagen University

% create matrix frame
A = zeros(N,N);

% normal distribution centered at d with std d/2
pd = makedist('Normal');

% START create frame (connect all vertices)
now = 1; 
NotVisited = 2:N;
while numel(NotVisited) > 1
    next = randsample(NotVisited,1,1);
    NotVisited = NotVisited(NotVisited ~= next);
    A(now,next) = 1;
    now = next;
end

A(now, NotVisited) = 1;
% END create frame

% number of edges
edges = N*d - 2*(N-1);

% START distribute edges among random nodes
idx = randi(N); % pick random node

while edges > 0
    entries = setdiff(find(A(idx,:)==0),idx); % find potential links of random node
    if ~isempty(entries) % if there are potential links...
        L = datasample(entries, min(min(max(1,round(abs(random(pd)))),round(edges/2)),numel(entries)),'Replace',false);  % ... pick a number of links based on distribution
        A(idx,L) = 1; % realize links
        edges = edges - 2*numel(L); % remove added links (edges) from the total number of edges required for mean degree d
    end
    
    idx = max(1,mod(idx+1,N+1));  % pick next node in line, if node > N pick instead node 1
    A = A + A'; % create symmetric matrix
    
end
% END distribute edges among random nodes

% finalize connection matrix. Make symmetric, normalize entries (to one)
% and remove diagonal values
A = A - diag(diag(A));
A = A + A';
A(A>1) = 1;


end

