function [ bp ] = createBoxplots( varargin )
% Lars Reiter Nielsen 
% CPH University - Science
% Creates multiple boxplots
% Input: multiple cells of (possible) dublicate elements

BPdata = varargin;
% if celldata-input, uncomment below
% BPdata = cellfun(@cell2mat,varargin,'UniformOutput',false);

C = cell2mat(BPdata);
 
 lengthArray = length(BPdata);
 lengthInternalArrays = cellfun('length',BPdata);
 
 grp = zeros(1,sum(lengthInternalArrays));
 for i = 1:lengthArray
     idx = find(grp == 0);
     grp(idx(1):idx(1)+lengthInternalArrays(i)-1) = i;    
 end
 
bp = boxplot(C,grp,'Whisker',4,'Symbol','ko','OutlierSize',1.5);
end

