classdef Galaxy < matlab.mixin.SetGet
    %GALAXY Summary of this class goes here
    
    properties
        stars Star
        
        starData
    end
    
    methods
        function obj = Galaxy(starData, optStarIds)
            obj.starData = starData;
            
            allStars = Star.empty(1,0);
            for(i=1:size(starData,1))
                row = starData(i,:);
                allStars(i) = Star(row(1), row(2), deg2rad(row(3)), deg2rad(row(4)), deg2rad(row(5)), deg2rad(row(6)), obj);
            end
            obj.stars = allStars;
            
            [~,Locb] = ismember(optStarIds,[obj.stars.id]);
            optStars = obj.stars(Locb);
            [optStars(:).isOptStar] = deal(true);
        end
        
        function IdStars = getStarsWithIDs(obj,ids)
            [~,Locb] = ismember(ids,[obj.stars.id]);
            IdStars = obj.stars(Locb);
        end
        
        function [validStars] = getValidColonizableStarsForAStar(obj, star, earliestColTime)
            validStars = findobj(obj.stars,'validColTarget',true, '-and', 'isColonized',false, '-and', 'isOptStar',true, '-and','-not', 'id',star.id);
%             validStars = findobj(obj.stars,'validColTarget',true, '-and', 'isColonized',false, '-and','-not', 'id',star.id);
%             bool = [validStars.inc] >= deg2rad(175);
%             validStars = validStars(bool);
                  
            % Find only stars in front of current star
%             [rVectCurStar, vVectCurStar] = getStarPositionKpcMyr(star.id, earliestColTime, obj.starData);
%             [rVectTgts, ~] = getStarPositionKpcMyr([validStars.id], earliestColTime, obj.starData);
%             
%             diffAtCurYear = rVectTgts - rVectCurStar;
%             vVectCurStar = repmat(vVectCurStar, 1, size(rVectTgts,2));
%             
%             dotProd = dot(vVectCurStar,diffAtCurYear,1);
%             normA = vecnorm(vVectCurStar,2,1);
%             normB = vecnorm(diffAtCurYear,2,1);
%             dAng = acosd(dotProd ./ (normA .* normB));
%             
%             validStars = validStars(dAng <= 90);
            
            % Whittle down base on distance at Myr 90
%             [rVectTgts, ~] = getStarPositionKpcMyr([validStars.id], 90, obj.starData);
%             
%             colStars = obj.getColonizedStars();
%             distMat = zeros(length(validStars), length(colStars));
%             for(i=1:length(colStars))
%                 rVectStar = getStarPositionKpcMyr(colStars(i).id, 90, obj.starData);
%                 dR = rVectTgts(1:2,:) - rVectStar(1:2,:);
%                 distMat(:,i) = rssq(dR,1);
%             end
%             
%             dMin = 0.5;
%             dMax = Inf;
%             
%             distances = min(distMat,[],2);
%             validStars = validStars(distances >= dMin & distances <= dMax);
            
%             rVectStar = getStarPositionKpcMyr(star.id, 90, obj.starData);
%             [rVectTgts, ~] = getStarPositionKpcMyr([validStars.id], 90, obj.starData);
%             dR = rVectTgts(1:2,:) - rVectStar(1:2,:);
%             distances = abs(rssq(dR,1) - 1);

%             [rVectCurStar, ~] = getStarPositionKpcMyr(star.id, earliestColTime, obj.starData);
%             [rVectTgts, ~] = getStarPositionKpcMyr([validStars.id], earliestColTime, obj.starData);
%             diffAtCurYear = rVectTgts - rVectCurStar;
%             normB = vecnorm(diffAtCurYear,2,1);
%             
%             numStars = 5000;
%             [~,I] = mink(normB,numStars);       
%             validStars = validStars(I);
        end
        
        function J = ipColObjFun(obj, x)
            x = unique(round(x),'stable');
            optStars = obj.getOptStars();
            optStars = optStars(x);
            
            exitColRun = false;
            
            colStars = obj.getColonizedStars();
            oldColStarIds = [colStars.id];
            curOptStarInd = 1;
            while(exitColRun == false)
                for(i=1:length(colStars))
                    for(j=1:3)
                        star1 = optStars(curOptStarInd);
                        colStars(i).colonizeFromThisStarToGivenStar(star1);

                        curOptStarInd = curOptStarInd + 1;
                    end
                end
                
                colStars = obj.getOptStarsThatHaventColonizedAnything();
                colStarIds = [colStars.id];
                
                Lia = ismember(oldColStarIds,colStarIds);
                if(all(Lia == 1))
                    exitColRun = true;
                end
            end
            
            J = obj.computeCurrentJ();
        end
                
        function colStars = getColonizedStars(obj)
            colStars = findobj(obj.stars, 'isColonized',true);
        end
        
        function optStars = getOptStars(obj)
            optStars = findobj(obj.stars, 'isOptStar',true);
        end
        
        function stars = getOptStarsThatHaventColonizedAnything(obj)
            optStars = findobj(obj.stars, 'isColonized',true, '-and', 'isOptStar',true);
            
            stars = Star.empty(1,0);
            for(i=1:length(optStars))
                star = optStars(i);
                subStars = obj.getStarsColonizedByStar(star);
                
                if(isempty(subStars))
                    stars(end+1) = star; %#ok<AGROW>
                end
            end
        end
        
        function stars = getStarsColonizedByStar(obj, star)
            stars = findobj(obj.stars, 'colonizedBy',star);
            
            if(isempty(stars))
                stars = Star.empty(1,0);
            end
        end
        
        function [fH,filename] = openSettlerShipOutputFile(obj)
            filename = sprintf('settlerShipARH_%s.txt',datestr(datetime(), 30));
            fH = fopen(filename,'w+');
        end
        
        function J = computeCurrentJ(obj)
            colStars = obj.getColonizedStars();
            colStarIds = [colStars.id];
            usedDv = [colStars.totalColDv];
            maxDv = [colStars.maxSettlerDv];
            [J,~,~] = computeMeritFunction(colStarIds, obj.starData, usedDv, maxDv);
        end
        
        function displayGraphOfColonization(obj)
            colStars = obj.getColonizedStars();
            
            srcsStars = [colStars.colonizedBy];
            
            srcs = [srcsStars.id];
            tgts = [colStars.id];
            
            nonColStarIds = setdiff(obj.starData(:,1),[srcs(:);tgts(:)]);
            nonColStarPos = getStarPositionKpcMyr(nonColStarIds,90, obj.starData);
            
            allNodes = unique([srcsStars(:); colStars(:)]);
            Rs = [allNodes.R];
            ThFs = [allNodes.thF];
            
            for(i=1:length(srcs))
                srcNames{i} = sprintf('Star %u',srcs(i)); %#ok<AGROW>
            end
            
            for(i=1:length(tgts))
                tgtNames{i} = sprintf('Star %u',tgts(i)); %#ok<AGROW>
            end
            
            weights = ones(size(srcs));
            
            allNames = cell(length(allNodes),1);
            for(i=1:length(allNodes))
                allNames{i} = sprintf('Star %u',allNodes(i).id);
            end
            [xdata,ydata] = pol2cart(ThFs,Rs);
            nodeTable = table(allNames,xdata(:),ydata(:),'VariableNames',{'Name','xdata','ydata'});
            
            G = digraph(srcNames,tgtNames,weights,nodeTable);
            
            hAx = axes(figure());
            hold on;
            plot(hAx,nonColStarPos(1,:),nonColStarPos(2,:),'.b');
            plot(hAx, G, '.r', 'XData',G.Nodes.xdata, 'YData',G.Nodes.ydata, 'MarkerSize',5);
            hold off;
            xlabel(hAx,'X [kpc]');
            ylabel(hAx,'Y [kpc]');
            axis(hAx,'equal');
            grid(hAx,'minor');
        end
    end
end