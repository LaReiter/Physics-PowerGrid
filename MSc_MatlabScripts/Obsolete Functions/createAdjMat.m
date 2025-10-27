function [ A ] = createAdjMat( N, p, fill)
% Creates 'pseudo-random' connection matrix with bin. distributed degree.
% N:    Number of nodes in the network
% p:    Probability parameter of binomial distribution. The network degree is
% directly proportional to p. Higher p makes for a more dense network, whereas
% lower p reduces the overall degree.
% fill:     If fill is set to 'true', the procedure will continue to add
% edges at each node drawn from the binomial distribution until each node
% has been visitied. If fill is set to 'false' the procedure will terminate
% as soon as every node is connected.
% NOTICE THAT THIS DOES NOT GUARANTEE AN OVERALL AVERAGE DEGREE N*p
% Code by Lars Reiter Nielsen @ Student at Copenhagen University

A = zeros(N,N); % initialize matrix 'frame'

MissingVisit = linspace(N,1,N);

% degree distribution. NOTE: Adjust second parameter to decrease/increase
% mean degree and degree variation
pd = makedist('Binomial',100,p); 

position = 1; % current node
MissingVisit = setdiff(MissingVisit,position); % current node is visited
oldposition = position; % last node
L = position; % initialize


if fill
    Visited = linspace(N,1,N); % vector of visited nodes
    Visited = setdiff(Visited,position); % current node is visited
end

while ~isempty(MissingVisit)

    while isempty(intersect(setdiff(L,[oldposition,position]),MissingVisit)) % continue until we reach a new node
        L = randperm(N, min(max(random(pd),1),N));
    end

%     A(position,union(L,oldposition))= 1; % make edges based on sample L
    A(position,L)=1;

    oldposition = position;
    
    if numel(intersect(setdiff(L,oldposition),MissingVisit)) == 1
        position = intersect(setdiff(L,oldposition),MissingVisit);
    else
        position = randsample(intersect(setdiff(L,oldposition),MissingVisit),1,1); % pick new position following random edge in L (avoid loop)
    end

    MissingVisit = setdiff(MissingVisit,L); % remove recently visited nodes
    
    if fill
        Visited = setdiff(Visited,position); % (position is visited) for later use **
    end

end

if fill
    for i=1:numel(Visited) %** add nodes to all non-visited nodes based on distribution
        A(Visited(i),randperm(N,min(max(random(pd),1),N))) = 1;
    end
end

A = A - diag(diag(A));
A = A + A' + tril(A)+tril(A)';
A(A>1) = 1;

        
% % %    CREATE MATRIC WITH AVERAGE DEGREE d
% % %         A = zeros(N,N); % define base matrix
% % %         while ~isConnected(A) && count < 5000
% % %             A = zeros(N,N); % define base matrix
% % %             Numbers = [ones(1,round(d*N/2)),zeros(1, N/2*(N-1)-round(d*N/2))];  % ones and zeroes
% % %             Numbers = Numbers(randperm(N/2*(N-1))); % permuted numbers with d*N ones
% % %             idx = 1;
% % %             for i=1:N-1
% % %                A(i,i+1:N) = Numbers(idx:idx+N-(1+i)); % add N-1 edges connecting every node
% % %               idx = idx + N-i;
% % %             end
% % %             A = A + A.'; % create symmetric matrix
% % %             count = count + 1;
% % %         end

end

