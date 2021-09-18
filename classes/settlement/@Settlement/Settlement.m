classdef Settlement < matlab.mixin.SetGet & matlab.mixin.CustomDisplay
    %Settlement Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fromNode(1,1) SettlementNode
        toNode(1,1) SettlementNode
        
        departTime(1,1) double
        flightTime(1,1) double
        
        context(1,1) SettlementContextEnum = SettlementContextEnum.SettlerShip
        
        totalDeltaV(1,1) double
        maxDeltaV(1,1) double
    end
    
    properties(Dependent)
        arriveTime(1,1) double
    end
    
    methods
        function obj = Settlement(fromNode, toNode, departTime, flightTime, context)
            obj.fromNode = fromNode;
            obj.toNode = toNode;
            obj.departTime = departTime;
            obj.flightTime = flightTime;
            obj.context = context;
        end
        
        function val = get.arriveTime(obj)
            val = obj.departTime + obj.flightTime;
        end
        
        function [deltaV1, deltaV2, deltaVTotal, trajOutsideBounds] = getDeltaVsForSettlement(obj, starData, maxDeltaV)
            starId1 = obj.fromNode.star.id;
            tStar1 = obj.departTime;
            starId2 = obj.toNode.star.id;
            tStar2 = obj.arriveTime;
            
            ids = [starId1; starId2];
            tMyr = [tStar1; tStar2];
            [~, starvVects] = getStarPositionKpcMyr(ids, tMyr, starData);
            
            %TM = 1
            [deltaV1(1), deltaV2(1), deltaVTotal(1)] = getDeltaVsForTM(obj, starData, starvVects, 1);
            
            %TM = -1
            [deltaV1(2), deltaV2(2), deltaVTotal(2)] = getDeltaVsForTM(obj, starData, starvVects, -1);
            
            if(deltaVTotal(1) < deltaVTotal(2))
                deltaV1 = deltaV1(1);
                deltaV2 = deltaV2(1);
                deltaVTotal = deltaVTotal(1);
            else
                deltaV1 = deltaV1(2);
                deltaV2 = deltaV2(2);
                deltaVTotal = deltaVTotal(2);
            end
            
            trajOutsideBounds = false;
            
            obj.totalDeltaV = deltaVTotal;
            obj.maxDeltaV = maxDeltaV;
        end
    end
    
    methods(Access=protected)
        function displayScalarObject(obj)
            fprintf('Settlement from Star %u to %u\n', obj.fromNode.star.id, obj.toNode.star.id);
        end
        
        function [deltaV1, deltaV2, deltaVTotal] = getDeltaVsForTM(obj, starData, starvVects, tm)
            starId1 = obj.fromNode.star.id;
            tStar1 = obj.departTime;
            starId2 = obj.toNode.star.id;
            tStar2 = obj.arriveTime;
            
            [~, vVect1, ~, vVect2] = gtocxStarKeplerLambert(starId1, tStar1, starId2, tStar2, tm, starData);
            
            deltaV1 = norm(vVect1 - starvVects(:,1));
            deltaV2 = norm(vVect2 - starvVects(:,2));
            deltaVTotal = deltaV1 + deltaV2;
        end
    end
end