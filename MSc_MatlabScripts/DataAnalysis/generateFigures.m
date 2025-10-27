
% Lars Reiter Nielsen
% CPH University

take1 = 0;
take2 = 0;
take3 = 0;
take4 = 0;
take5 = 0;
take6 = 1;

%% CBPGBPfunc (DR)
while take1
    figure()
for i=1:22
   CBPGBPfunc(DynamicalArray{1,1}{1,3,i},DynamicalArray{1,1}{1,1,i},0);
   disp(i);
end

for i=1:21
   CBPGBPfunc(DynamicalArray{1,2}{1,3,i},DynamicalArray{1,2}{1,1,i},0);
   disp(i);
end

for i=1:20
   CBPGBPfunc(DynamicalArray2{1,3}{1,3,i},DynamicalArray2{1,3}{1,1,i},0);
   disp(i);
end

for i=1:49
   CBPGBPfunc(DynamicalArray3{1,4}{1,3,i},DynamicalArray3{1,4}{1,1,i},0);
   disp(i);
end

for i=1:49
   CBPGBPfunc(DynamicalArray4{1,5}{1,3,i},DynamicalArray4{1,5}{1,1,i},0);
   disp(i);
end

for i=1:49
   CBPGBPfunc(DynamicalArray4{1,6}{1,3,i},DynamicalArray4{1,6}{1,1,i},0);
   disp(i);
end

for i=1:49
   CBPGBPfunc(DynamicalArray4{1,7}{1,3,i},DynamicalArray4{1,7}{1,1,i},0);
   disp(i);
end

for i=1:49
   CBPGBPfunc(DynamicalArray4{1,8}{1,3,i},DynamicalArray4{1,8}{1,1,i},0);
   disp(i);
end

for i=1:49
   CBPGBPfunc(DynamicalArray4{1,9}{1,3,i},DynamicalArray4{1,9}{1,1,i},0);
   disp(i);
end

for i=1:49
   CBPGBPfunc(DynamicalArray5{1,10}{1,3,i},DynamicalArray5{1,10}{1,1,i},0);
   disp(i);
end

for i=1:49
   CBPGBPfunc(DynamicalArray5{1,11}{1,3,i},DynamicalArray5{1,11}{1,1,i},0);
   disp(i);
end

for i=1:49
   CBPGBPfunc(DynamicalArray5{1,12}{1,3,i},DynamicalArray5{1,12}{1,1,i},0);
   disp(i);
end

take1=0;
    
end


%% CBPGBPfunc (TR)
while take2
    for j=2:12
        for i=1:3
            CBPGBPfunc(TransArray{1,j}{1,3,i},TransArray{1,j}{1,1,i},1);
        end
    end
    for j=1:12
        CBPGBPfunc(TransArray2{1,j}{1,3,1},TransArray2{1,j}{1,1,1}, 1);
    end
    for i=1:3
            CBPGBPfunc(TransArray3{1,1}{1,3,i},TransArray3{1,1}{1,1,i}, 1);
    end
    for j=2:4
        for i=1:3
            CBPGBPfunc(TransArray4{1,j}{1,3,i},TransArray4{1,j}{1,1,i}, 1);
        end
    end
    take2=0;
end

%% Boxplots (DR)
while take3
    BP1 = [];
    for i=1:22
        BP1 = cat(2,BP1,DynamicalArray{1,1}{1,1,i});
    end
    createBoxplots(BP1, DynamicalArray{2,2}{2}, DynamicalArray2{2,3}{2}, DynamicalArray3{2,4}{2},DynamicalArray4{2,5}{2},DynamicalArray4{2,6}{2},DynamicalArray4{2,7}{2},DynamicalArray4{2,8}{2},DynamicalArray4{2,9}{2},DynamicalArray5{2,10}{2},DynamicalArray5{2,11}{2},DynamicalArray5{2,12}{2})
    take3=0;
end

%% Normalplots (DR)
while take4
    
    MagnitudeMean1 = zeros(1,22);
    for i = 1:22
        MagnitudeMean1(i) = mean(DynamicalArray{1,1}{1,1,i});
    end
    
    MagnitudeMean2 = zeros(1,21);
    for i=1:21
       MagnitudeMean2(i) = mean(DynamicalArray{1,2}{1,1,i});
    end

    MagnitudeMean3 = zeros(1,20);
    for i=1:20
       MagnitudeMean3(i) = mean(DynamicalArray2{1,3}{1,1,i});
    end

    MagnitudeMean4 = zeros(1,49);
    for i=1:49
       MagnitudeMean4(i) = mean(DynamicalArray3{1,4}{1,1,i});
    end

    MagnitudeMean5 = zeros(1,49);
    for i=1:49
        MagnitudeMean5(i) = mean(DynamicalArray4{1,5}{1,1,i});
    end

    MagnitudeMean6 = zeros(1,49);
    for i=1:49
        MagnitudeMean6(i) = mean(DynamicalArray4{1,6}{1,1,i});
    end

    MagnitudeMean7 = zeros(1,49);
    for i=1:49
        MagnitudeMean7(i) = mean(DynamicalArray4{1,7}{1,1,i});
    end

    MagnitudeMean8 = zeros(1,49);
    for i=1:49
        MagnitudeMean8(i) = mean(DynamicalArray4{1,8}{1,1,i});
    end

    MagnitudeMean9 = zeros(1,49);
    for i=1:49
       MagnitudeMean9(i) = mean(DynamicalArray4{1,9}{1,1,i});
    end

    MagnitudeMean10 = zeros(1,49);
    for i=1:49
       MagnitudeMean10(i) = mean(DynamicalArray5{1,10}{1,1,i});
    end

    MagnitudeMean11 = zeros(1,49);
    for i=1:49
        MagnitudeMean11(i) = mean(DynamicalArray5{1,11}{1,1,i});
    end

    MagnitudeMean12 = zeros(1,49);
    for i=1:49
       MagnitudeMean12(i) = mean(DynamicalArray5{1,12}{1,1,i});
    end
    take4=0;
    MagnitudesArray = {MagnitudeMean1,MagnitudeMean2,MagnitudeMean3,MagnitudeMean4,MagnitudeMean5,MagnitudeMean6,MagnitudeMean7,MagnitudeMean8,MagnitudeMean9,MagnitudeMean10,MagnitudeMean11,MagnitudeMean12};
    for j=1:3:12
        createNormalplots(MagnitudesArray{j},MagnitudesArray{j+1},MagnitudesArray{j+2});
        figure()
    end
end

%% Boxplots (TR)
while take5
    BP1=[];
    BP1 = cat(2,BP1,TransArray2{1,1}{1,1,1});
    BP1 = cat(2,BP1, TransArray3{1,1}{1,1,1});
    BP1 = cat(2,BP1, TransArray3{1,1}{1,1,2});
    BP1 = cat(2,BP1, TransArray3{1,1}{1,1,3});
    BP1 = cat(2,BP1, TransArray3{1,1}{1,1,4});
    
    BP2 = [];
    BP2 = cat(2,BP2,TransArray{1,2}{1,1,1});
    BP2 = cat(2,BP2,TransArray{1,2}{1,1,2});
    BP2 = cat(2,BP2,TransArray{1,2}{1,1,3});
    BP2 = cat(2, BP2, TransArray2{1,2}{1,1,1});
    BP2 = cat(2, BP2, TransArray4{1,2}{1,1,1});
    BP2 = cat(2, BP2, TransArray4{1,2}{1,1,2});
    BP2 = cat(2, BP2, TransArray4{1,2}{1,1,3});
    
    BP3 = [];
    BP4 = [];
    BP5 = [];
    BP6 = [];
    BP7 = [];
    BP8 = [];
    BP9 = [];
    BP10 = [];
    BP11 = [];
    BP12 = [];
    for i = 1:3
        BP12 = cat(2,BP12, TransArray4{1,12}{1,1,i});
        BP11 = cat(2,BP11, TransArray4{1,11}{1,1,i});
        BP10 = cat(2,BP10, TransArray4{1,10}{1,1,i});
        BP9 = cat(2,BP9, TransArray4{1,9}{1,1,i});
        BP8 = cat(2,BP8, TransArray4{1,8}{1,1,i});
        BP7 = cat(2,BP7, TransArray4{1,7}{1,1,i});
        BP6 = cat(2,BP6, TransArray4{1,6}{1,1,i});
        BP5 = cat(2,BP5, TransArray4{1,5}{1,1,i});
        BP4 = cat(2,BP4,TransArray{1,4}{1,1,i});
        BP4 = cat(2,BP4,TransArray4{1,4}{1,1,i});
        BP3 = cat(2,BP3,TransArray{1,3}{1,1,i});
        BP3 = cat(2, BP3, TransArray4{1,3}{1,1,i});
    end
    BP3 = cat(2, BP3, TransArray2{1,3}{1,1,1});
    BP4 = cat(2,BP4,TransArray2{1,4}{1,1,1});
    BP4 = cat(2,BP4,TransArray2{1,4}{1,1,2});
    BP12 = cat(2,BP12, TransArray2{1,12}{1,1,1});
    BP12 = cat(2,BP12, TransArray2{1,12}{1,1,2});
    BP11 = cat(2,BP11, TransArray2{1,11}{1,1,1});
    BP10 = cat(2,BP10, TransArray2{1,10}{1,1,1});
    BP10 = cat(2,BP10, TransArray2{1,10}{1,1,2});
    BP9 = cat(2,BP9, TransArray2{1,9}{1,1,1});
    BP8 = cat(2,BP8, TransArray2{1,8}{1,1,1});
    BP7 = cat(2,BP7, TransArray2{1,7}{1,1,1});
    BP6 = cat(2,BP6, TransArray2{1,6}{1,1,1});
    BP5 = cat(2,BP5, TransArray2{1,5}{1,1,1});
    take5=0;
    
    createBoxplots(BP1,BP2,BP3,BP4,BP5,BP6,BP7,BP8,BP9,BP10,BP11,BP12)
end

%% Normalplots (TR)
while take6
    BP1=[];
    BP1 = cat(2,BP1,mean(TransArray2{1,1}{1,1,1}));
    BP1 = cat(2,BP1, mean(TransArray3{1,1}{1,1,1}));
    BP1 = cat(2,BP1, mean(TransArray3{1,1}{1,1,2}));
    BP1 = cat(2,BP1, mean(TransArray3{1,1}{1,1,3}));
    BP1 = cat(2,BP1, mean(TransArray3{1,1}{1,1,4}));
    
    BP2 = [];
    BP2 = cat(2,BP2,mean(TransArray{1,2}{1,1,1}));
    BP2 = cat(2,BP2,mean(TransArray{1,2}{1,1,2}));
    BP2 = cat(2,BP2,mean(TransArray{1,2}{1,1,3}));
    BP2 = cat(2, BP2, mean(TransArray2{1,2}{1,1,1}));
    BP2 = cat(2, BP2, mean(TransArray4{1,2}{1,1,1}));
    BP2 = cat(2, BP2, mean(TransArray4{1,2}{1,1,2}));
    BP2 = cat(2, BP2, mean(TransArray4{1,2}{1,1,3}));
    
    BP3 = [];
    BP4 = [];
    BP5 = [];
    BP6 = [];
    BP7 = [];
    BP8 = [];
    BP9 = [];
    BP10 = [];
    BP11 = [];
    BP12 = [];
    for i = 1:3
        BP12 = cat(2,BP12, mean(TransArray4{1,12}{1,1,i}));
        BP11 = cat(2,BP11, mean(TransArray4{1,11}{1,1,i}));
        BP10 = cat(2,BP10, mean(TransArray4{1,10}{1,1,i}));
        BP9 = cat(2,BP9, mean(TransArray4{1,9}{1,1,i}));
        BP8 = cat(2,BP8, mean(TransArray4{1,8}{1,1,i}));
        BP7 = cat(2,BP7, mean(TransArray4{1,7}{1,1,i}));
        BP6 = cat(2,BP6, mean(TransArray4{1,6}{1,1,i}));
        BP5 = cat(2,BP5, mean(TransArray4{1,5}{1,1,i}));
        BP4 = cat(2,BP4,mean(TransArray{1,4}{1,1,i}));
        BP4 = cat(2,BP4,mean(TransArray4{1,4}{1,1,i}));
        BP3 = cat(2,BP3,mean(TransArray{1,3}{1,1,i}));
        BP3 = cat(2, BP3, mean(TransArray4{1,3}{1,1,i}));
    end
    BP3 = cat(2, BP3, mean(TransArray2{1,3}{1,1,1}));
    BP4 = cat(2,BP4,mean(TransArray2{1,4}{1,1,1}));
    BP4 = cat(2,BP4,mean(TransArray2{1,4}{1,1,2}));
    BP12 = cat(2,BP12, mean(TransArray2{1,12}{1,1,1}));
    BP12 = cat(2,BP12, mean(TransArray2{1,12}{1,1,2}));
    BP11 = cat(2,BP11, mean(TransArray2{1,11}{1,1,1}));
    BP10 = cat(2,BP10, mean(TransArray2{1,10}{1,1,1}));
    BP10 = cat(2,BP10, mean(TransArray2{1,10}{1,1,2}));
    BP9 = cat(2,BP9, mean(TransArray2{1,9}{1,1,1}));
    BP8 = cat(2,BP8, mean(TransArray2{1,8}{1,1,1}));
    BP7 = cat(2,BP7, mean(TransArray2{1,7}{1,1,1}));
    BP6 = cat(2,BP6, mean(TransArray2{1,6}{1,1,1}));
    BP5 = cat(2,BP5, mean(TransArray2{1,5}{1,1,1}));
    take6=0;
    
    BasinStabArray = {BP1,BP2,BP3,BP4,BP5,BP6,BP7,BP8,BP9,BP10,BP11,BP12};
    for j=1:3:12
        createNormalplots(BasinStabArray{j},BasinStabArray{j+1},BasinStabArray{j+2});   
    end
end