clc; clear variables; format long g; close all;

addpath(genpath('data'));
addpath(genpath('methods'));
addpath(genpath('solutions'));

% starData = readmatrix('data/stars.txt');
starData = dlmread('data/stars.txt',',',1,0);
filteredSd = starData(starData(:,3) >= 179.9 & starData(:,2) <= 15,:);
% filteredSd = starData;

% x0 = [0, deg2rad(90), deg2rad(0), 0.1, ...
%       1, deg2rad(90), deg2rad(0), 0.1];
% x0 = [      5.220378
%             1.695140323614
%        -0.0403466732051036
%                0.012008367
%                44.27924488
%           4.54687730717959
%        0.00333567320510375
%                0.028667206];
ub = [45, 2*pi, pi/2, 0.2, ...
      45, 2*pi, pi/2, 0.2];
lb = [0, 0, -pi/2, 0.01, ...
      1, 0, -pi/2, 0.01];

% fVect = objFun(x0, [], filteredSd);

fun = @(x, xdata) objFun(x, xdata, filteredSd);
% options = optimoptions(@lsqcurvefit, 'Display','iter', 'MaxFunctionEvaluations',Inf, 'MaxIterations',5000);
% [x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(fun,x0,[],zeros(size(fVect)),lb,ub,options);

%some initial guesses
nvars = length(ub);
pop = (ub-lb).*rand(1000,length(ub)) + lb;

% Create OPTI Object
funNorm = @(x) objFunNorm(x, [], filteredSd);
opts = optiset('solver','nomad', 'display','iter', 'maxiter',10000, 'maxfeval',1E6, 'maxtime',3600);
Opt = opti('fun',funNorm,'bounds',lb,ub,'options',opts);

%parallel run
numRuns = 7;
xx = nan([nvars,numRuns]);
fval = nan([1,numRuns]);
exitflag = nan([1,numRuns]);
info = cell(1,numRuns);
% parfor(i=1:numRuns)
%     [xx(:,i),fval(i),exitflag(i),info{i}] = solve(Opt,pop(i,:));
% end

%single run
[xx,fval,exitflag,info] = solve(Opt,pop(1,:));

function fVect = objFun(x, ~, starData)
    tBurn1 = x(1);
    azBurn1 = x(2);
    elBurn1 = x(3);
    magBurn1 = x(4);
    
    dtBurn1Burn2 = x(5);
    azBurn2 = x(6);
    elBurn2 = x(7);
    magBurn2 = x(8);
    
    %Step size
    step = 0.001;
    
    %get init position @ Sol
    [rVectKpc, vVectKpcMyr] = getStarPositionKpcMyr(0, tBurn1, starData);

    %apply DV 1
    [dvx,dvy,dvz] = sph2cart(azBurn1,elBurn1,magBurn1);
    vVectKpcMyr = vVectKpcMyr + [dvx,dvy,dvz]';
    
    %propagate to burn 2
    t0 = tBurn1;
    dt = dtBurn1Burn2;
    [t1,y1, ~,~,~] = propagateBody(t0,rVectKpc,vVectKpcMyr,dt, step, []);

    %apply DV 2
    rVectKpc = y1(end,1:3)';
    vVectKpcMyr = y1(end,4:6)';
    [dvx,dvy,dvz] = sph2cart(azBurn2,elBurn2,magBurn2);
    vVectKpcMyr = vVectKpcMyr + [dvx,dvy,dvz]';
    
    %propagate to 90 Myr
    t0 = tBurn1 + dtBurn1Burn2;
    dt = 90 - t0;
    [t2,y2, ~,~,~] = propagateBody(t0,rVectKpc,vVectKpcMyr,dt, step, []);
    
    %concat arrays
    t = vertcat(t1,t2);
    y = vertcat(y1,y2);
    
    %compute star positions at all times for all stars (w/o Sol #0)
    starData(starData(:,1) == 0, :) = [];
    starIds = starData(:,1);
    A = allcomb(t,starIds);
    
    [rVectStars, ~] = getStarPositionKpcMyr(A(:,2), A(:,1), starData);
    
    %compute residuals all all time steps
    rVects = y(:,1:3)'; 
    a = repmat(rVects, length(starIds), 1);
    rVects = reshape(a(:),3,numel(a)/3);
    
    res = rVects - rVectStars;
    
	fVect = rssq(res,1)';
    fVect = reshape(fVect,length(t),length(starIds));
    fVect = min(fVect,[],1)';
end

function f = objFunNorm(x, ~, starData)
    fVect = objFun(x, [], starData);
%     fVect(fVect > 1) = 1;
%     
%     if(isempty(fVect))
%         f = 10;
%     else
%         f = norm(fVect);
%     end
    
    f = norm(fVect);
end