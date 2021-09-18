clc; clear variables; clear format long g; close all;

addpath(genpath('data'));
addpath(genpath('methods'));
addpath(genpath('solutions'));
addpath(genpath('studies'));

starData = dlmread('data/stars.txt',',',1,0);
optStarIds = unique(csvread('solutions/mothership_team/optStars/nomadoptstarssolfile_350n_J332.0.txt'));

% optStarIds = optStarIds(1:150);
times = 0:1:90;
nodes = combvec(optStarIds',times)';
numNodes = size(nodes,1);

[rVectStars, vVectStars] = getStarPositionKpcMyr(nodes(:,1), nodes(:,2), starData);

ppm = ParforProgMon('Scanning...', numNodes);

starIds = nodes(:,1);
times = nodes(:,2);
srcs = nan([numNodes,numNodes]);
tgts = nan([numNodes,numNodes]);
weights = nan([numNodes,numNodes]);
nodenames = cell(numNodes,1);

profile on;
for(i=1:5) %numNodes
    disp(i);
    starId0 = starIds(i);
    departTime = times(i);
    departRVect = rVectStars(:,i);
    departVVect = vVectStars(:,i);
    
    targetNodes = nodes(nodes(:,2) >= (departTime+1) & nodes(:,1) ~= starId0, :); 
    depNodeId = i;
    
    srcsAllJ = nan(1,numNodes);
    tgtsAllJ = nan(1,numNodes);
    weightsAllJ = nan(1,numNodes);
    for(j=1:length(targetNodes))
        starId1 = targetNodes(j,1);
        arriveTime = targetNodes(j,2);
        
        arrNodeId = find(nodes(:,1) == starId1 & nodes(:,2) == arriveTime,1);
        arriveRVect = rVectStars(:,arrNodeId);
        arriveVVect = vVectStars(:,arrNodeId);
        
        tof = arriveTime-departTime;
        x = [departTime, tof];
        [dvTot, dv1, dv2, ~, ~] = dvForOneShip(x, departRVect, departVVect, arriveRVect, arriveVVect, false);
        
        if(dv1 < 175 && dv2 < 175 && dvTot < 400)
            srcsAllJ(j) = depNodeId;
            tgtsAllJ(j) = arrNodeId;
            weightsAllJ(j) = tof;
            
        end
    end
    
    srcs(i,:) = srcsAllJ;
    tgts(i,:) = tgtsAllJ;
    weights(i,:) = weightsAllJ;
    
    nodenames{i} = sprintf('Star %u @ T = %0.3f Myr', starId0,departTime);

%     ppm.increment(); %#ok<PFBNS>
end
profile viewer;

% srcs = srcs(:);
% tgts = tgts(:);
% weights = weights(:);
% 
% srcs(isnan(srcs)) = [];
% tgts(isnan(tgts)) = [];
% weights(isnan(weights)) = [];
% 
% G = digraph(srcs,tgts,weights,nodenames);
% plot(G,'XData',nodes(:,1),'YData',nodes(:,2));
% xlim([0,max(nodes(:,1))]);
% ylim([0 90]);
% xlabel('Star ID');
% ylabel('Departure Time from Star [MYr]');
% grid on;


        function [dvTot, dv1, dv2, dv1Vect, dv2Vect] = dvForOneShip(x, rVect0, vVect0, rVect1, vVect1, isHighFidelity)
            t0 = x(1);
            dt = x(2);
            t1 = t0+dt;

            [dvTot1, dv11, dv12, dv11Vect, dv12Vect] = computeLambertArcDv(rVect0, vVect0, t0, rVect1, vVect1, t1, 1, isHighFidelity);
            [dvTot2, dv21, dv22, dv21Vect, dv22Vect] = computeLambertArcDv(rVect0, vVect0, t0, rVect1, vVect1, t1, -1, isHighFidelity);

            if(dvTot1 < dvTot2)
                dvTot = dvTot1; 
                dv1 = dv11;
                dv2 = dv12;
                dv1Vect = dv11Vect;
                dv2Vect = dv12Vect;
            else
                dvTot = dvTot2;
                dv1 = dv21;
                dv2 = dv22;
                dv1Vect = dv21Vect;
                dv2Vect = dv22Vect;
            end
        end

        function [dvTot, dv1, dv2, dv1Vect, dv2Vect] = computeLambertArcDv(rVect0, vVect0, t0, rVect1, vVect1, t1, tm, isHighFidelity)
%             [~, vVect0] = getStarPositionKpcMyr(starId0, t0, starData);
%             [~, vVect1] = getStarPositionKpcMyr(starId1, t1, starData);

            if(isHighFidelity)
%                 [~, vVect0t, ~, vVect1t] = gtocxStarLambert(starId0, t0, starId1, t1, tm, starData);
                [vVect0t, vVect1t] = gtocxRVectLambert(rVect0, t0, rVect1, t1, tm);
            else
%                 [~, vVect0t, ~, vVect1t] = gtocxStarKeplerLambert(starId0, t0, starId1, t1, tm, starData);
                [vVect0t, vVect1t] = gtocxRVectKeplerLambert(rVect0, t0, rVect1, t1, tm);
            end

            dv1Vect = (vVect0t-vVect0)*(1/kmS2KpcMyr);
            dv1 = norm(dv1Vect);
            
            dv2Vect = (vVect1-vVect1t)*(1/kmS2KpcMyr);
            dv2 = norm(dv2Vect);
            
            dvTot = (dv1 + dv2);
        end