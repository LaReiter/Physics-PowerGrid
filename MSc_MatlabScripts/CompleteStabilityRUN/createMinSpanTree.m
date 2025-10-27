function [ A ] = createMinSpanTree( N, g )
%createAdjMat creates a N x N connection matrix with N-1 edges (MST)
%   Input parameters:
%   N:  Dimension of the square (symmetric) connection/adjacency matrix
%   (N = number of nodes) Mininum span tree has N-1 edges
%   g: node granularity, number between 1 and N
% Code by Lars Reiter Nielsen @ Student at Copenhagen University

% create matrix frame
A = zeros(N,N);

% Visited nodes have entry 1 otherwise 0
Visited = 2:N;

% current node
here = 1;

while ~isempty(Visited)
    sample = datasample(Visited, min(randi(length(Visited)),round(N/g)), 'Replace', false);
%     sample = datasample(Visited,randi(length(Visited)), 'Replace', false);
    A(here, sample) = 1;
    Visited = setdiff(Visited, sample);
    here = datasample(sample, 1);
end

A = A + A';
A(A>1) = 1;


end

