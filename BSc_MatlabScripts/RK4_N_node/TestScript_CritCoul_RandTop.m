% Test script
clc;

% initial conditions
N0 = 5;
NG = 2;
edges = N0;
alpha = 0.1;
dt = 1e-1;
time = 200;


pmaxes = [];
rep = 0;
repMax = 5;

while rep < repMax;
    % matr = random_graph_modified(N0,1,edges);
    w0 = zeros(N0,1);
    t0 = zeros(N0,1);
    p0 = -ones(N0,1);
    gens = datasample(1:N0,NG,'Replace',false);
    p0(gens) = (N0-NG)/NG;

    [p, Rinf] = findKc(time, dt, alpha, mat, p0, w0, t0)
    
    if Rinf>0.9
        pmaxes(end+1)=p;
        rep = rep+1;
    end
end
disp(Rinf)
disp(p)