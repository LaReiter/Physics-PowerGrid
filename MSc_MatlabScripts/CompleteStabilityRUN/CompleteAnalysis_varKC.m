clc; clear all;
% complete stability analysis
% Lars Reiter Nielsen
% DynamicalArray{1,i} contains DATAarray of type i (check below)
% DynamicalArray{2,i} contains CELLarray of type i
% - DATAarray: DATAarray{1,j,g} - g is the current grid (g=1..G), j=1 is
% the MagnitudeList of grid g, j = 2 is the connection matrix of grid g, j
% =3 is the power vector of grid g, j = 4 is the frequency of the
% magnitudes of magnitudelist, i.e. frequency(i) = MagnitudeList(i)/N_g
% where N_g is the number of nodes of grid g
%
% - CELLarray: CELLarray{1,i} - i = 1 is a magnitude vector
% (1+advance:advance:magMax), i=2 is an occurence vector MagFreq which contains all measured
% magnitudes among all grids of the grid-type and includes potential dublicates, that
% is,MagFreq(1:length(DATAarray{1,1,1}) for instance includes measured
% magnitudes of grid g=1 for all nodes (in order i=1...N_g)
% all simulated grids, i = 3 is a meanMagVec consisting of the magnitude
% means of each simulated grid (used to make normal distribution plot)
TransArray = cell(1,12);
DynamicalArray = cell(2,12);
t = 1/4*1e-3;
KC = 4;

% DYNAMICAL

% type 1

% try
%    [A,B] = DynamicTestFunc_ST_VarKC(1, 'C', t);
%    DynamicalArray{1,1} = A;
%    DynamicalArray{2,1} = B;
% catch
%    disp(strcat('Derror ',strcat(num2str(1),'C')))
% end

try
   [A2,B2] = DynamicTestFunc_ST_VarKC(1, 'U', t);
   DynamicalArray{1,2} = A;
   DynamicalArray{2,2} = B;
catch
   disp(strcat('Derror ',strcat(num2str(1),'U')))
end

try
   [A,B] = DynamicTestFunc_ST_VarKC(1, 'R', t);
   DynamicalArray{1,3} = A;
   DynamicalArray{2,3} = B;
catch
   disp(strcat('Derror ',strcat(num2str(1),'R')))
end


% type 2

try
   [A,B] = DynamicTestFunc_ST_VarKC(2, 'C', t);
   DynamicalArray{1,4} = A;
   DynamicalArray{2,4} = B;
catch
   disp(strcat('Derror ',strcat(num2str(2),'C')))
end

try
   [A,B] = DynamicTestFunc_ST_VarKC(2, 'U', t);
   DynamicalArray{1,5} = A;
   DynamicalArray{2,5} = B;
catch
   disp(strcat('Derror ',strcat(num2str(2),'U')))
end

try
   [A,B] = DynamicTestFunc_ST_VarKC(2, 'R', t);
   DynamicalArray{1,6} = A;
   DynamicalArray{2,6} = B;
catch
   disp(strcat('Derror ',strcat(num2str(2),'R')))
end

% type 3

try
   [A,B] = DynamicTestFunc_ST_VarKC(3, 'C', t);
   DynamicalArray{1,7} = A;
   DynamicalArray{2,7} = B;
catch
   disp(strcat('Derror ',strcat(num2str(3),'C')))
end

try
   [A,B] = DynamicTestFunc_ST_VarKC(3, 'U', t);
   DynamicalArray{1,8} = A;
   DynamicalArray{2,8} = B;
catch
   disp(strcat('Derror ',strcat(num2str(3),'U')))
end

try
   [A,B] = DynamicTestFunc_ST_VarKC(3, 'R', t);
   DynamicalArray{1,9} = A;
   DynamicalArray{2,9} = B;
catch
   disp(strcat('Derror ',strcat(num2str(3),'R')))
end

% type 4

try
   [A,B] = DynamicTestFunc_ST_VarKC(4, 'C', t);
   DynamicalArray{1,10} = A;
   DynamicalArray{2,10} = B;
catch
   disp(strcat('Derror ',strcat(num2str(4),'C')))
end

try
   [A,B] = DynamicTestFunc_ST_VarKC(4, 'U', t);
   DynamicalArray{1,11} = A;
   DynamicalArray{2,11} = B;
catch
   disp(strcat('Derror ',strcat(num2str(4),'U')))
end

try
   [A,B] = DynamicTestFunc_ST_VarKC(4, 'R', t);
   DynamicalArray{1,12} = A;
   DynamicalArray{2,12} = B;
catch
   disp(strcat('Derror ',strcat(num2str(4),'R')))
end





% TRANSIENT

% % type 1
% 
% % try
% %    A = TransientTestFunc(1,'C',t,KC);
% %    TransArray{1,1} = A;
% % catch
% %    disp(strcat('Terror ',strcat(num2str(1),'C')))
% % end
% % t = 1e-3;
% try
%    A = TransientTestFunc(1, 'U', t,KC);
%    TransArray{1,2} = A;
% catch
%    disp(strcat('Terror ',strcat(num2str(1),'U')))
% end
% 
% try
%    A = TransientTestFunc(1, 'R', t,KC);
%    TransArray{1,3} = A;
% catch
%    disp(strcat('Terror ',strcat(num2str(1),'R')))
% end
% 
% % type 2
% 
% try
%    A = TransientTestFunc(2, 'C', t,KC);
%    TransArray{1,4} = A;
% catch
%    disp(strcat('Terror ',strcat(num2str(2),'C')))
% end
% 
% try
%    A = TransientTestFunc(2, 'U', t,KC);
%    TransArray{1,5} = A;
% catch
%    disp(strcat('Terror ',strcat(num2str(2),'U')))
% end
% 
% try
%    A = TransientTestFunc(2, 'R', t,KC);
%    TransArray{1,6} = A;
% catch
%    disp(strcat('Terror ',strcat(num2str(2),'R')))
% end
% 
% % type 3
% 
% try
%    A = TransientTestFunc(3, 'C', t,KC);
%    TransArray{1,7} = A;
% catch
%    disp(strcat('Terror ',strcat(num2str(3),'C')))
% end
% 
% try
%    A = TransientTestFunc(3, 'U', t,KC);
%    TransArray{1,8} = A;
% catch
%    disp(strcat('Terror ',strcat(num2str(3),'U')))
% end
% 
% try
%    A = TransientTestFunc(3, 'R', t,KC);
%    TransArray{1,9} = A;
% catch
%    disp(strcat('Terror ',strcat(num2str(3),'R')))
% end
% 
% % type 4
% 
% try
%    A = TransientTestFunc(4, 'C', t,KC);
%    TransArray{1,10} = A;
% catch
%    disp(strcat('Terror ',strcat(num2str(4),'C')))
% end
% 
% try
%    A = TransientTestFunc(4, 'U', t,KC);
%    TransArray{1,11} = A;
% catch
%    disp(strcat('Terror ',strcat(num2str(4),'U')))
% end
% 
% try
%    A = TransientTestFunc(4, 'R', t,KC);
%    TransArray{1,12} = A;
% catch
%    disp(strcat('Terror ',strcat(num2str(4),'R')))
% end
   