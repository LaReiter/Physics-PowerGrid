function [ np ] = createNormalplots( varargin )
% Lars Reiter Nielsen 
% CPH University - Science
% Creates multiple normalplots
% Input: multiple arrays of measurements

meanData = varargin;
% if inputdata is celldata, uncomment code snippet below
% meanMagData = cellfun(@cell2mat,varargin,'UniformOutput',false);

% set maximum depending on data (one-node basin stability or magnitudes)
x = 0:0.01:1.5;
maxy=1;
for i =1:length(varargin)
   norm = normpdf(x, mean(meanData{i}),std(meanData{i}));
   maxy = max(maxy,max(norm));
end

figure()
hold on
patch([1 1.5 1.5 1], [0 0 maxy*1.1 maxy*1.1], [0.9 0.9 0.9])
for i =1:length(varargin)
   norm = normpdf(x, mean(meanData{i}),std(meanData{i}));
   plot(x, norm,'LineWidth',2);
end
legend({'Unaccesable', 'Complete','Uniform','Random'});
ylim([0,maxy*1.1])
% choose legend according to data
% 
% meanMagMat = cell2mat(meanMagData')'; % arrange magnitudes in lists
% 
% meanMag = mean(meanMagMat);
% stdMag = std(meanMagMat);
end

