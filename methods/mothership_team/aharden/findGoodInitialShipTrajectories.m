clc; clear variables; format long g; close all;

addpath(genpath('data'));
addpath(genpath('methods'));
addpath(genpath('solutions'));
% parfevalOnAll(gcp('nocreate'), @warning, 0, 'off');

% starData = readmatrix('data/stars.txt');
starData = dlmread('data/stars.txt',',',1,0);

% optStarIds = unique(csvread('data/solN11208.csv'));
% filteredSd = starData(optStarIds+1,:);
% filteredSd = [starData(1,:);filteredSd];
filteredSd = starData(starData(:,3)>179,:);

popSize = 10000;
numSettlePods = 2;

lbTmpl = [ 1,  1,  2];
ubTmpl = [55, 55,  size(filteredSd,1)-1];

fun = @(x) initShipTrajObjFcnCoarse(x, filteredSd);
nvars = 3*3*numSettlePods;
lb = repmat(lbTmpl,[1,nvars/3]);
lb(1:3*numSettlePods:nvars) = 0;

ub = repmat(ubTmpl,[1,nvars/3]);
ub(1:3*numSettlePods:nvars) = 10;

IntCon = [3:3:nvars];
nonlcon = @(x) initShipTrajNonlconCoarse(x, filteredSd);

A = zeros(3,nvars);
A(1,setdiff(0*numSettlePods+1:1:3*numSettlePods, IntCon)) = 1;
A(2,setdiff(3*numSettlePods+1:1:6*numSettlePods, IntCon)) = 1;
A(3,setdiff(6*numSettlePods+1:1:9*numSettlePods, IntCon)) = 1;

b = zeros(3,1);
b(1:3,:) = 50;

pop = (ub-lb).*rand(popSize,nvars) + lb;
pop(:,IntCon) = round(pop(:,IntCon));
nonlconTest = nonlcon(pop(1,:));

% options = optimoptions(@ga, 'Display','iter', 'UseParallel',false, 'PopulationSize',size(pop,1), 'PlotFcn',{'gaplotscores','gaplotbestf'}, 'InitialPopulationMatrix',pop, 'MaxGenerations',100, 'MaxStallGenerations',20);
% [xx,fval,exitflag,output,population,scores] = ga(fun,nvars,A,b,[],[],lb,ub,nonlcon,IntCon,options);

% Options
% xtype = repmat('C',[1,nvars]);
% xtype(IntCon) = 'I';
% opts = optiset('solver','nomad', 'display','iter', 'maxiter',10000, 'maxfeval',1E6, 'maxtime',3600);
% 
% % Create OPTI Object
% nlcon = nonlcon;
% nlrhs = zeros(size(nonlconTest));
% nle = -1*ones(size(nonlconTest));
% %'nlmix',nlcon,nlrhs,nle,
% Opt = opti('fun',fun,'ineq',A,b,'bounds',lb,ub,'xtype',xtype,'options',opts);

% warning('off','all');
% options = optimoptions(@patternsearch,'ScaleMesh',true, 'Display','iter', 'UseParallel',false, 'SearchFcn',{@searchlhs,5,15}, 'InitialMeshSize',1E-5, 'Cache','on', 'CacheSize',1E7, 'UseCompletePoll',true, 'UseCompleteSearch',true, 'MaxFunctionEvaluations',Inf, 'MaxIterations',Inf, 'MaxTime',3600);
% [xx,fval] = patternsearch(fun,pop(100,:),A,b,[],[],lb,ub,nonlcon,options);
% warning('on','all');

% Solve the MINLP problems
numRuns = 32*10;
xx = nan([numRuns,nvars]);
fval = nan([1,numRuns]);
exitflag = nan([1,numRuns]);
info = cell(1,numRuns);
nonlconRun = nan(numRuns,length(nonlconTest));
options = optimoptions(@patternsearch,'ScaleMesh',false, 'Display','iter', 'UseParallel',false, 'PollMethod','MADSPositiveBasis2N', 'SearchFcn',{@searchlhs,5,15}, 'InitialMeshSize',1000, 'Cache','on', 'CacheSize',1E7, 'UseCompletePoll',false, 'UseCompleteSearch',false, 'MaxFunctionEvaluations',Inf, 'MaxIterations',Inf, 'MaxTime',3600);
warning('off','all');
parfor(i=1:numRuns)
%     [xx(:,i),fval(i),exitflag(i),info{i}] = solve(Opt,pop(i,:));
%     nonlconRun{i} = nonlcon(xx(:,i));
    [xx(i,:),fval(i),exitflag(i),info{i}] = patternsearch(fun,pop(i,:),A,b,[],[],lb,ub,nonlcon,options);
    nonlconRun(i,:) = nonlcon(xx(i,:))';
end
warning('on','all');

% [xx,fval,exitflag,info] = solve(Opt,pop(1,:));
% disp(max(nonlcon(xx)))

bool = max(nonlconRun,[],2) <= 1E-6;

if(any(bool))
    xx = xx(bool,:);
    fval = fval(bool);

    [~,I] = min(fval);
else
    [~,I]=min(max(nonlconRun,[],2));
end

x = xx(I,:);

[~, msCsvMatrix] = fun(x);

starIdList = round(x(IntCon));
starIdList = filteredSd(starIdList+1,1);
dts = x(setdiff(1:nvars,IntCon));
save('mothershipRuns.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% funF = @(x) initShipTrajObjFcnFine(x, starIdList, starData);
% x0F = dts;
% 
% AF = zeros(3,length(x0F));
% AF(1,0*numSettlePods+1:1:2*numSettlePods) = 1;
% AF(2,2*numSettlePods+1:1:4*numSettlePods) = 1;
% AF(3,4*numSettlePods+1:1:6*numSettlePods) = 1;
% 
% bF = zeros(3,1);
% bF = b;
% 
% lbF = 1*ones(1,length(x0F));
% lbF(1:numSettlePods*2:length(x0F)) = 0;
% 
% ubF = 90*ones(1,length(x0F));
% ubF(1:numSettlePods*2:length(x0F)) = 90;
% 
% nonlconF = @(x) initShipTrajNonlconFine(x, starIdList, starData);
% nonlconTestF = nonlconF(x0F);
% 
% warning('off','all');
% options = optimoptions(@fmincon, 'Display','iter', 'UseParallel',true, 'ConstraintTolerance',1E-6);
% [xF,fvalF] = fmincon(funF,x0F,AF,bF,[],[],lbF,ubF,nonlconF,options);
% warning('on','all');
% 
% % warning('off','all');
% % options = optimoptions(@patternsearch, 'Display','iter', 'UseParallel',true, 'ConstraintTolerance',1E-6);
% % [xF,fvalF] = patternsearch(funF,x0F,AF,bF,[],[],lbF,ubF,nonlconF,options);
% % warning('on','all');
% 
% [~, msCsvMatrix] = funF(xF);

% % nlcon = nonlcon;
% % nlrhs = zeros(size(nonlconTest));
% % nle = -1*ones(size(nonlconTest));
% % 
% % %Options
% % opts = optiset('solver','ipopt', 'display','iter');
% % 
% % % Solve Problem
% % Opt = opti('fun',fun,'ineq',A,b,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub, 'options',opts);
% % [x2,fval2] = solve(Opt,x0);