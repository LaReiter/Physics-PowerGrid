 function [ P ] = createPowerVec( A, ConFrac, MaxPow, Uniform)
% createPowerVec creates vector respresenting consumed and generated power.
% Input parameters:
%   A: Connection matrix (N x N)
%   ConFrac: Fraction of consumers (between 0 and 1). 
%   MaxPow: Upper bound for the maximum nodewise consumed/generated power.
%   Uniform: If 'true', the power is distributed evenly among consumers and
%   evenly among generators
% Output:  power vector P (N x 1) representing consumers (P<0) and generators
% (P>0).
% Code by Lars Reiter Nielsen @ Student at Copenhagen University

N = length(A);
P = zeros(N,1);  % initialize N x 1 vector

noCon = ceil(ConFrac*N);
noGen = N-noCon;

Consumers = randperm(N, noCon); % random consumers
Generators = setdiff(linspace(N,1,N),Consumers); % random generators

    if ConFrac <= 0.5
        if ~Uniform
            P(Consumers) = -MaxPow*rand(1,noCon);
            randNums = rand(1,noGen);
            P(Generators) = sum(-P(Consumers))*(1/sum(randNums)*randNums);
        else
            P(Consumers) = -MaxPow*repmat(rand(1,1),1,noCon);
            P(Generators) = sum(-P(Consumers))*repmat(1/noGen,1,noGen);
        end
    else
        if ~Uniform
            P(Generators) = MaxPow*rand(1,noGen);
            randNums = rand(1,noCon);
            P(Consumers) = sum(-P(Generators))*(1/sum(randNums)*randNums);
        else
            P(Generators) = MaxPow*repmat(rand(1,1),1,noGen);
            P(Consumers) = sum(-P(Generators))*repmat(1/noCon,1,noCon);
        end        
    end

end

