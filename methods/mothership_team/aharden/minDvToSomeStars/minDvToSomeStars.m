clc; clear variables; clear format long g; close all; fclose('all');

addpath(genpath('.'));

starData = dlmread('data/stars.txt',',',1,0);
optStarIds = unique(csvread('data/solN11208.csv'));
filteredSd = starData(optStarIds+1,:);
filteredSd = filteredSd(filteredSd(:,3)>179,:);

numStars = 2;
numRuns = 1;

funLf = @(x) objFun(x, false, filteredSd, starData);
nonlconLf = @(x) nonlconFun(x, funLf);

funHf = @(x) objFun(x, false, filteredSd, starData);
nonlconHf = @(x) nonlconFun(x, funHf);

lbTmpl = [1,  1];
ubTmpl = [90, size(filteredSd,1)];

lb = repmat(lbTmpl,1,numStars);
ub = repmat(ubTmpl,1,numStars);

pop = (ub-lb).*rand(numRuns,length(lb)) + lb;
pop(:,[2:2:length(lb)]) = round(pop(:,[2:2:length(lb)]));

A = zeros(1,length(lb));
A(:,1:2:length(lb)) = 1;
b = [90];

% options = optimoptions(@patternsearch,'ScaleMesh',false, 'Display','iter', 'UseParallel',false, 'PollMethod','MADSPositiveBasis2N', 'SearchFcn',{@searchlhs,5,15}, 'InitialMeshSize',0.1, 'Cache','on', 'CacheSize',1E7, 'UseCompletePoll',true, 'UseCompleteSearch',true, 'MaxFunctionEvaluations',Inf, 'MaxIterations',Inf, 'MaxTime',3600);
% [xLf, fvalLf] = patternsearch(funLf,pop(1,:),A,b,[],[],lb,ub,nonlconLf,options);

options = optimoptions(@fmincon, 'MaxIterations',Inf, 'MaxFunctionEvaluations',Inf);
problemLf = createOptimProblem('fmincon','objective',funLf, 'x0',pop(1,:), 'lb',lb, 'ub',ub, 'Aineq',A, 'bineq',b, 'nonlcon',nonlconLf, 'options',options);
gs = GlobalSearch('Display','iter');
[xLf, fvalLf] = gs.run(problemLf);
[~, eachDvLf] = funLf(xLf);

% lb(2:2:length(lb)) = xLf(2:2:length(lb));
% ub(2:2:length(ub)) = xLf(2:2:length(ub));
% 
% problemHf = createOptimProblem('fmincon','objective',funHf, 'x0',xLf, 'lb',lb, 'ub',ub, 'Aineq',A, 'bineq',b);
% % [xHf, fvalHf] = gs.run(problemHf);
% [xHf, fvalHf] = fmincon(problemHf);
% [~, eachDvHf] = funHf(xHf);

function [c, ceq] = nonlconFun(x, objFun)
    [~, eachDv1] = objFun(x);
    c = eachDv1(2:end) - (0.01 - 1E-6);
    
    ceq = [];
end

function [f, eachDv1, eachDv2] = objFun(x, isHighFidelity, filteredSd, starData)
    t0 = 0;
    starId0 = 0;
    
    [rVect0, vVect0] = getStarPositionKpcMyr(starId0, t0, starData);
    
    totalDv = 0;
    eachDv1 = [];
    eachDv2 = [];
    for(i=1:length(x)/2)
        dt = x(2*i-1);
        
        starInd = round(x(2*i));
        starId1 = filteredSd(starInd,1);
        
        [dvTot, dv1, dv2, dv1Vect, dv2Vect, rVect1t, vVect1t, trajOutsideBounds] = dvForOneShip(t0, rVect0, vVect0, dt, starId1, isHighFidelity, starData);
        totalDv = totalDv + dv1;
        eachDv1(i) = dv1; %#ok<AGROW>
        eachDv2(i) = dv2; %#ok<AGROW>
        
        t0 = t0 + dt;
        rVect0 = rVect1t;
        vVect0 = vVect1t;
    end
    
    f = totalDv;
end

        function [dvTot, dv1, dv2, dv1Vect, dv2Vect, rVect1t, vVect1t, trajOutsideBounds] = dvForOneShip(t0, rVect0, vVect0, dt, starId1, isHighFidelity, starData)
            t1 = t0+dt;

            [dvTot1, dv11, dv12, dv11Vect, dv12Vect, rVect11t, vVect11t, trajOutsideBounds1] = computeLambertArcDv(rVect0, vVect0, t0, starId1, t1, 1, isHighFidelity, starData);
            [dvTot2, dv21, dv22, dv21Vect, dv22Vect, rVect12t, vVect12t, trajOutsideBounds2] = computeLambertArcDv(rVect0, vVect0, t0, starId1, t1, -1, isHighFidelity, starData);

            if(dvTot1 < dvTot2)
                dvTot = dvTot1; 
                dv1 = dv11;
                dv2 = dv12;
                dv1Vect = dv11Vect;
                dv2Vect = dv12Vect;
                rVect1t =  rVect11t;
                vVect1t = vVect11t;
                trajOutsideBounds = trajOutsideBounds1;
            else
                dvTot = dvTot2;
                dv1 = dv21;
                dv2 = dv22;
                dv1Vect = dv21Vect;
                dv2Vect = dv22Vect;
                rVect1t =  rVect12t;
                vVect1t = vVect12t;
                trajOutsideBounds = trajOutsideBounds2;
            end
        end

        function [dvTot, dv1, dv2, dv1Vect, dv2Vect, rVect1, vVect1t, trajOutsideBounds] = computeLambertArcDv(rVect0, vVect0, t0, starId1, t1, tm, isHighFidelity, starData)
%             [~, vVect0] = getStarPositionKpcMyr(starId0, t0, starData);
            [rVect1, vVect1] = getStarPositionKpcMyr(starId1, t1, starData);

            if(isHighFidelity)
%                 [~, vVect0t, ~, vVect1t, trajOutsideBounds] = gtocxStarLambert(starId0, t0, starId1, t1, tm, starData);
                [vVect0t, vVect1t, trajOutsideBounds] = gtocxRVectLambert(rVect0, t0, rVect1, t1, tm);
            else
%                 [~, vVect0t, ~, vVect1t] = gtocxStarKeplerLambert(starId0, t0, starId1, t1, tm, starData);
                [vVect0t, vVect1t] = gtocxRVectKeplerLambert(rVect0, t0, rVect1, t1, tm);
                trajOutsideBounds = false;
            end

            dv1Vect = (vVect0t-vVect0)*(1/kmS2KpcMyr);
            dv1 = norm(dv1Vect);
            
            dv2Vect = (vVect1-vVect1t)*(1/kmS2KpcMyr);
            dv2 = norm(dv2Vect);
            
            dvTot = (dv1 + dv2);
        end