classdef SettlementGraph < matlab.mixin.SetGet
    %SettlementGraph Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nodes(:,1) SettlementNode
        settlements(:,1) Settlement
        
        G
        
        starData StarsData
    end
    
    methods
        function obj = SettlementGraph(initNode, starData)
            obj.starData = starData;
            obj.addNode(initNode);
        end
        
        function addNode(obj, node)
            obj.nodes(end+1) = node;
        end
        
        function addSettlement(obj, settlement)
            obj.settlements(end+1) = settlement;
            settlement.fromNode.addDepartingSettlement(settlement);
        end
        
        function removeSettlement(obj, settlement)
            obj.settlements(obj.settlements == settlement) = [];
            settlement.fromNode.removeDepartingSettlement(settlement);
        end
        
        function [earliestAvailableNode, settlementContext] = getEarliestAvailableNode(obj)           
            availableNodes = obj.nodes([obj.nodes.remainingSettlements] > 0 & [obj.nodes.settledOn] < 87);
            [~,I] = min([availableNodes.settledOn]);
            earliestAvailableNode = availableNodes(I);
            
            if(isempty(earliestAvailableNode))
                earliestAvailableNode = SettlementNode.empty(1,0);
                settlementContext = SettlementContextEnum.SettlerShip;
                
                return;
            end
            
%             for(i=1:length(obj.nodes))
%                 node = obj.nodes(i);
%                 subSettlements = obj.settlements([obj.settlements.fromNode] == node);
%                 numSettlementsFromThisNode = numel(subSettlements);
%                 
%                 if(node.star.id == 0)
%                     N = 5;
%                 else
%                     N = 3;
%                 end
%                 
%                 if(numSettlementsFromThisNode < N) %need to handle fast ships and motherships at some point
%                     if(isempty(earliestAvailableNode))
%                         earliestAvailableNode = node;
%                     else
%                         if(earliestAvailableNode.settledOn > node.settledOn)
%                             earliestAvailableNode = node;
%                         end
%                     end
%                 end
%             end
            
            if(not(isempty(earliestAvailableNode)))
                if(earliestAvailableNode.star.id == 0)
                    subSettlements = obj.settlements([obj.settlements.fromNode] == earliestAvailableNode);
                    numSettlementsFromThisNode = numel(subSettlements);

                    if(numSettlementsFromThisNode < 2)
                        settlementContext = SettlementContextEnum.FastShip;
                    else
                        settlementContext = SettlementContextEnum.Mothership;
                    end
                else
                    settlementContext = SettlementContextEnum.SettlerShip;
                end
            end
        end
        
        function obsInfo = getGalaxyStarObservationInfo(obj)
            obsInfo = obj.starData.getTemplateObsMatrix();
            
            settledStarIds = obj.getSettledStarIds();
            
            for(i=1:length(settledStarIds))
                settledStarId = settledStarIds(i);
                obsInfo(obsInfo(:,1) == settledStarId, 4) = 1;
            end
        end
        
        function node = getNodeForDepartureStar(obj, star)
            fromNodes = [obj.nodes];
            node = fromNodes([fromNodes.star] == star);
        end
        
        function settledStarIds = getSettledStarIds(obj)
            settledStars = [obj.nodes.star];
            settledStarIds = [settledStars.id];
        end
        
        function J = getObjFunValue(obj)
            settledStarIds = obj.getSettledStarIds();
            
            if(numel(settledStarIds) > 1 && numel(obj.settlements) > 0)
                actualDvs = NaN(numel(obj.settlements),1);
                maxDvs = NaN(numel(obj.settlements),1);
                for(i=1:length(obj.settlements))
                    actualDvs(i) = obj.settlements(i).totalDeltaV;
                    maxDvs(i) = obj.settlements(i).maxDeltaV;
                end

                [J,~,~] = computeMeritFunction(settledStarIds, obj.starData.starsData, actualDvs(:), maxDvs(:));
            else
                J = 0;
            end
        end
        
        function plotSettlementGraph(obj)
            hAx = axes(figure());
            hold(hAx,'on');
            axis(hAx,'equal');
            axis(hAx,'tight');
            grid(hAx,'minor');
            xlim(hAx,[-32,32]);
            ylim(hAx,[-32,32]);
            
            for(i=1:length(obj.settlements))
                settlement = obj.settlements(i);
                star1Id = settlement.fromNode.star.id;
                tStar1 = settlement.fromNode.settledOn;
                
                star2Id = settlement.toNode.star.id;
                tStar2 = settlement.toNode.settledOn;
                
                [starrVects, ~] = getStarPositionKpcMyr([star1Id;star2Id], [tStar1;tStar2], obj.starData.starsData);
                plot3(hAx, starrVects(1,:), starrVects(2,:), starrVects(3,:), 'ro');
                
                posDiff = starrVects(:,2) - starrVects(:,1);
                quiver3(hAx, starrVects(1,1), starrVects(2,1), starrVects(3,1), posDiff(1), posDiff(2), posDiff(3), 'r', 'AutoScale','off');
            end
        end
    end
end