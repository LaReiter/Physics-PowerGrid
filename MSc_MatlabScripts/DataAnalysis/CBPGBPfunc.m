function [ ] = CBPGBPfunc( Pvec, Measurements, Trans)
% Lars Reiter Nielsen
% Given Pvec determines the value of CBP and GBP and from this calculate
% the skewness:
% skew(CBP) = (mean(CBP) - median(CBP))/std(CBP) aswell as for GBP
% makes a single point in a plot where "Measurements" is a list of either 
% magnitudes or single-node basin stabilities. The size of the point
% indicates the average magnitude strength OR the average basin stability
% of the grid in question. If Trans = true, the pointsize is scaled to
% basin one-node stability
% program displays a plot with skewness(delta) = x, and skewness(mu) = y
% where each point has size proportional to measurement strength

Pc = [1/5,1/2, 1, 2, 6]; % power of consumers
Pg = [1/5,33/20, 31/10, 91/20, 6]; % power of generators

consumers = sort(abs(Pvec(Pvec < 0))); % sort consumer power elements (ascending)
generators = sort(Pvec(Pvec> 0)); % sort generator power elements (asc.)

conElem = unique(consumers); %  vector of unique consumer power elemennts
genElem = unique(generators); % vector of unique generator power elements

c =histc(consumers,conElem);  % vector of number OF unique consumer power elements
g =histc(generators,genElem); % vector of number OF unique generator power elements

% delta = 1/25*g.*genElem/sum(1/25*g.*genElem);
% mu = 1/25*c.*conElem/sum(1/25*c.*conElem);

% determine CBP (mu)
mu = zeros(1,5);
[dummyvar, indicesConsumers] = intersect(Pc,consumers); % find size of consumer elements and correpsonding indices
mu(indicesConsumers') = c/sum(c); % rebuild CBP


% determine GBP (delta)
delta= zeros(1,5);
indicesGenerators = zeros(1,length(genElem));  % indices for GBP
for i = length(genElem):-1:1
    [dummyvar2, idxGen] = min(abs(Pg-genElem(i))); % the lowest difference between gen power element and Pg vector elements must yield proper index
    indicesGenerators(i) = idxGen;
    Pg(idxGen) = 1e10;
end
delta(indicesGenerators) = g/sum(g); % rebuild GBP

dummyGenVar1 = zeros(1,5);
dummyGenVar2 = zeros(1,5);
dummyConVar1 = zeros(1,5);
dummyConVar2 = zeros(1,5);
for i=1:length(indicesGenerators)
    dummyGenVar1(indicesGenerators(i)) = genElem(i);
    dummyGenVar2(indicesGenerators(i)) = g(i);
end
for i=1:length(indicesConsumers)
   dummyConVar1(indicesConsumers(i)) = conElem(i);
   dummyConVar2(indicesConsumers(i)) = c(i);
end


delta = 1/25*dummyGenVar2.*dummyGenVar1/sum(1/25*dummyGenVar2.*dummyGenVar1);
mu = 1/25*dummyConVar2.*dummyConVar1/sum(1/25*dummyConVar2.*dummyConVar1);
    

dummyDelta = [];
dummyMu = [];
for i = 1:5
   dummyDelta = cat(2,dummyDelta, repmat(i,1,round(delta(i)*100)));
   dummyMu = cat(2, dummyMu, repmat(i,1,round(mu(i)*100)));
end
% orderVec = 1:5; % tiny to huge
% skewdelta = (mean(delta.*orderVec) - median(delta.*orderVec))/std(delta.*orderVec);
% skewmu = (mean(mu.*orderVec)-median(mu.*orderVec))/std(mu.*orderVec);
skewdelta = (mean(dummyDelta) - median(dummyDelta))/std(dummyDelta);
skewmu = (mean(dummyMu)-median(dummyMu))/std(dummyMu);

if Trans
    mmax = 1; % one-node basin max
    col = 'r';
    markmax = 8;
else
    mmax = 8; % magnitude max (semi-arbitrary)
    col='b';
    markmax = 12;
end
markersize = (mean(Measurements)/mmax);
plot(skewdelta,skewmu,'o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',markersize*markmax)

hold on
end

