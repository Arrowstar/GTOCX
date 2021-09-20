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
        function obj = Settlement(fromNode, toNode, context)
            obj.fromNode = fromNode;
            obj.toNode = toNode;
            obj.context = context;
        end
        
        function val = get.arriveTime(obj)
            val = obj.departTime + obj.flightTime;
        end
        
%         function [deltaV1, deltaV2, deltaVTotal, trajOutsideBounds] = getDeltaVsForSettlement(obj, starData, maxDeltaV)
%             starId1 = obj.fromNode.star.id;
%             tStar1 = obj.departTime;
%             starId2 = obj.toNode.star.id;
%             tStar2 = obj.arriveTime;
%             
%             ids = [starId1; starId2];
%             tMyr = [tStar1; tStar2];
%             [~, starvVects] = getStarPositionKpcMyr(ids, tMyr, starData);
%             
%             %TM = 1
%             [deltaV1(1), deltaV2(1), deltaVTotal(1)] = getDeltaVsForTM(obj, starData, starvVects, 1);
%             
%             %TM = -1
%             [deltaV1(2), deltaV2(2), deltaVTotal(2)] = getDeltaVsForTM(obj, starData, starvVects, -1);
%             
%             if(deltaVTotal(1) < deltaVTotal(2))
%                 deltaV1 = deltaV1(1);
%                 deltaV2 = deltaV2(1);
%                 deltaVTotal = deltaVTotal(1);
%             else
%                 deltaV1 = deltaV1(2);
%                 deltaV2 = deltaV2(2);
%                 deltaVTotal = deltaVTotal(2);
%             end
%             
%             trajOutsideBounds = false;
%             
%             obj.totalDeltaV = deltaVTotal;
%             obj.maxDeltaV = maxDeltaV;
%         end
        
        function [deltaV1, deltaV2, deltaVTotal, trajOutsideBounds, exitflag] = getOptimizedDeltaV(obj, starData, settledOnTime)
            fun = @(x) obj.deltaVObjFun(x, settledOnTime, starData);
            x0 = [3.0001, 5];
            lb = [3, 0.01];
            ub = [10, 40];
            A = [1 1];
            b = 90 - settledOnTime;
            nonlcon = @(x) obj.deltaVObjNonlcon(x, settledOnTime, starData);
            options = optimoptions(@fmincon, 'Display','none', 'MaxIterations',200);
            
%             [x,~,exitflag,~] = fmincon(fun,x0,[],[],[],[],lb,ub,nonlcon,options);
            problem = createOptimProblem('fmincon', 'objective',fun, 'x0',x0, 'Aineq',A, 'bineq',b, 'lb',lb, 'ub',ub, 'nonlcon',nonlcon, 'options',options);
            ms = MultiStart();
            [x,~,exitflag,~] = run(ms, problem,5);
            
            [deltaV1, deltaV2, deltaVTotal] = obj.getDeltaVsForX(x, settledOnTime, starData);
            trajOutsideBounds = false;
            
            waitTime = x(1);
            tof = x(2);
            
            tStar1 = settledOnTime + waitTime;
            tStar2 = tStar1 + tof;
            
            obj.departTime = tStar1;
            obj.flightTime = tof;
            
            obj.totalDeltaV = deltaVTotal;
            obj.maxDeltaV = obj.context.maxTotalDv;
            obj.toNode.settledOn = tStar2;
        end
    end
    
    methods(Access = private)
        function f = deltaVObjFun(obj, x, settledOnTime, starData)
            waitTime = x(1);
            tof = x(2);
            
            [~, ~, deltaVTotal] = obj.getDeltaVsForX(x, settledOnTime, starData);
            f = deltaVTotal + (waitTime + tof)/10;
        end
        
        function [c, ceq] = deltaVObjNonlcon(obj, x, settledOnTime, starData)
            [deltaV1, deltaV2, deltaVTotal] = obj.getDeltaVsForX(x, settledOnTime, starData);
            
            c(1) = deltaV1 - obj.context.maxIndivDv;
            c(2) = deltaV2 - obj.context.maxIndivDv;
            c(3) = deltaVTotal - obj.context.maxTotalDv;
            
            ceq = [];
        end
        
        function [deltaV1, deltaV2, deltaVTotal] = getDeltaVsForX(obj, x, settledOnTime, starData)
            starId1 = obj.fromNode.star.id;
            starId2 = obj.toNode.star.id;
            
            waitTime = x(1);
            tof = x(2);
            
            tStar1 = settledOnTime + waitTime;
            tStar2 = tStar1 + tof;
            
            [~, starvVects] = getStarPositionKpcMyr([starId1;starId2], [tStar1;tStar2], starData);
            
            %Delta-v TM = 1
            [~, vVect1, ~, vVect2] = gtocxStarKeplerLambert(starId1, tStar1, starId2, tStar2, 1, starData);
            
            deltaV1(1) = norm(vVect1 - starvVects(:,1));
            deltaV2(1) = norm(vVect2 - starvVects(:,2));
            deltaVTotal(1) = deltaV1(1) + deltaV2(1);
            
            %Delta-v TM = -1
            [~, vVect1, ~, vVect2] = gtocxStarKeplerLambert(starId1, tStar1, starId2, tStar2, -1, starData);
            deltaV1(2) = norm(vVect1 - starvVects(:,1));
            deltaV2(2) = norm(vVect2 - starvVects(:,2));
            deltaVTotal(2) = deltaV1(2) + deltaV2(2);
            
            if(deltaVTotal(1) < deltaVTotal(2))
                deltaV1 = deltaV1(1);
                deltaV2 = deltaV2(1);
                deltaVTotal = deltaVTotal(1);
            else
                deltaV1 = deltaV1(2);
                deltaV2 = deltaV2(2);
                deltaVTotal = deltaVTotal(2);
            end
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