function [ A ] = createBiConGraph( N, E )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nodes = [1];
E0 = E;
E1 = E0-1;

while ~isempty(nodes) || E0 > E1
    A = zeros(N,N);

    for i = 1:N
        for j = i+1:N
            A(i,j) = randi([0,1],1);
            A(j,i) = A(i,j);
        end
    end
    
     G = graph(A);
     
     E1 = sum(sum(triu(A)));

    [edges,nodes] = biconncomp(G);
end



E0 = E1 - E0;

% reached = 0;

while E0 > 0
   
   degrees = sum(A);
   
   Node2List = [];
   Node1List = find(degrees>2);
   
   while isempty(Node2List)
       node1 = datasample(Node1List,1);

       Node2List = intersect(Node1List, find(A(node1,:)>0));
   end
   node2 = datasample(Node2List,1);
       
   A(node1,node2)=0;
   A(node2,node1)=0;

   E0 = E0-1;

end

