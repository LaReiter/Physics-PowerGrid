function [ bool ] = isConnected( A )
% isConnected checks if there is a path between all nodes in A.
% A:    Symmetric square matrix A representing a nodal network. If i,j are nodes
% then A(i,j) > 0 is interpreted as a connection and A(i,j) = 0 means no connection.
% Code by Lars Reiter Nielsen @ Student at Copenhagen University


g = digraph(A);
bins = conncomp(g,'Type','weak');
bool = all(bins==1);


% % % Initialization
% % dim = length(A); % matrix dimension
% % Visited = zeros(1,dim); % initially none of the nodes are visited
% % 
% % % all nodes are connected if and only if node 1 is connected to all nodes
% % % We check that node 1 is connected to the rest:
% % if sum(sum(A)) == 0
% %     bool = 0;
% % else
% %     for i=1:dim
% %         B = A^i;    
% %         Visited(find(B(1,:) >0)) = 1; % check if i-th degree of seperation
% %     end
% %     
% %         % if all nodes are visited
% %     if sum(Visited) == dim
% %         bool = 1;
% %     else
% %         bool = 0;
% % end
% % end



end

