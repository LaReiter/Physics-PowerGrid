 function [ ] = dispGraph( A, P)
% Given connection matrix A (N x N) creates and displays a undirected graph 
% based on the adjacency matrix A
% Code by Lars Reiter Nielsen @ Student at Copenhagen University

h = plot(graph(A));
generators = find(P>0);
highlight(h, generators, 'NodeColor','r', 'Marker', 's')

end

