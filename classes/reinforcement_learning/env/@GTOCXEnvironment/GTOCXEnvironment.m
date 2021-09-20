classdef GTOCXEnvironment < rl.env.MATLABEnvironment
    %GTOCXEnvironment: Template for defining custom environment in MATLAB.    
    
    %% Properties (set properties' attributes accordingly)
    properties
        % Specify and initialize environment's necessary properties    
        starsData StarsData
        settlementGraph SettlementGraph
        
        %current state parameters
        currentTime(1,1) double = 0;
        currentNode(1,1) SettlementNode
        currentSettlementContext(1,1) SettlementContextEnum = SettlementContextEnum.FastShip
    end
    
    properties(Access = protected)
        % Initialize internal flag to indicate episode termination
        IsDone(1,1) logical = false        
    end

    %% Necessary Methods
    methods              
        % Contructor method creates an instance of the environment
        % Change class name and constructor name accordingly
        function obj = GTOCXEnvironment(starsDataLoc)
            starsData = StarsData(starsDataLoc);
            
            solStar = starsData.getStarForId(0);           
            initNode = SettlementNode(solStar,5);
            initNode.settledOn = -3;
            settlementGraph = SettlementGraph(initNode, starsData);
            
            % Initialize Observation settings
            ObservationInfo = GTOCXEnvironment.getInitObjInfo(settlementGraph);
            
            % Initialize Action settings   
            ActionInfo = GTOCXEnvironment.getInitActionInfo(settlementGraph);
            
            % The following line implements built-in functions of RL env
            obj = obj@rl.env.MATLABEnvironment(ObservationInfo,ActionInfo);
            
            %set current state properties
            obj.currentNode = initNode;
            obj.currentTime = 0;
            obj.currentSettlementContext = SettlementContextEnum.FastShip;
            
            %set env properties
            obj.starsData = starsData;
            obj.settlementGraph = settlementGraph;
        end
        
        % Apply system dynamics and simulates the environment with the 
        % given action for one step.
        function [Observation,Reward,IsDone,LoggedSignals] = step(obj,Action)
            LoggedSignals = [];
            
            actInfo = obj.getActionInfo();
            
            bool = Action <= actInfo.LowerLimit;
            Action(bool) = actInfo.LowerLimit(bool);
%             
            bool = Action >= actInfo.UpperLimit;
            Action(bool) = actInfo.UpperLimit(bool);
            
            %Process Star ID
%             a = obj.settlementGraph.getGalaxyStarObservationInfo();
%             allIds = a(:,1);
            
%             starId = round(Action(1));
%             if(starId > max(allIds))
%                 starId = max(allIds);
%             elseif(starId < min(allIds))
%                 starId = min(allIds);
%             end

            starR = Action(1);
            starThetaF = Action(2);
            starIds = obj.starsData.getStarClosestTo(starR, starThetaF);

            %%% Handle colonization
            settledStarIds = obj.settlementGraph.getSettledStarIds();
            I = find(~ismember(starIds, settledStarIds),1,'first');
            starId = starIds(I);
            if(~ismember(starId, settledStarIds))

                %Get Arrive Star
                arriveStar = obj.starsData.getStarForId(starId);

                %Get depart and arrive nodes
                departNode = obj.currentNode;
                arriveNode = SettlementNode(arriveStar, 3);

                %create settlement
                context = obj.currentSettlementContext;
                settlement = Settlement(departNode, arriveNode, context);

                %Evaluate DVs                
                starData = obj.starsData.starsData;
                [~, ~, ~, ~, exitflag] = settlement.getOptimizedDeltaV(starData, obj.currentTime);
                                
                if(exitflag > 0)
                    penalty = 0.1;
                    
                    %Add settlement and settlement node to graph if delta-v
                    %is satisfied
                    obj.settlementGraph.addNode(arriveNode);
                    obj.settlementGraph.addSettlement(settlement);
                else
                    penalty = -0.1;
                end
                
                [earliestAvailableNode, newSettlementContext] = obj.settlementGraph.getEarliestAvailableNode();

                if(not(isempty(earliestAvailableNode)))
                    time = earliestAvailableNode.settledOn;
                    star = earliestAvailableNode.star;
                    Observation = obj.getObsInfoForTimeAndStarId(time, star.id);

%                     disp(time);
                    
                    obj.currentTime = time;
                    obj.currentNode = earliestAvailableNode;
                    obj.currentSettlementContext = newSettlementContext;
                    
                    if(time >= 87)
                        IsDone = 1;
                    else
                        IsDone = 0;
                    end

                    Reward = penalty;
                    if(IsDone)
                        Reward = Reward + obj.settlementGraph.getObjFunValue();
                    end
                else
                    IsDone = 1;
                    
                    Reward = -penalty;
                    Reward = Reward + obj.settlementGraph.getObjFunValue();
                end
            else
                Observation = obj.getObsInfoForTimeAndStarId(obj.currentTime, obj.currentNode.star.id);
                Reward = obj.settlementGraph.getObjFunValue();
                IsDone = 1;
            end
            
            % (optional) use notifyEnvUpdated to signal that the 
            % environment has been updated (e.g. to update visualization)
            notifyEnvUpdated(obj);
        end
        
        % Reset environment to initial state and output initial observation
        function InitialObservation = reset(obj)
            initTime = 0;
            initStarId = 0;
            
            solStar = obj.starsData.getStarForId(initStarId);
            initNode = SettlementNode(solStar, 5);
            initNode.settledOn = -3;
            obj.settlementGraph = SettlementGraph(initNode, obj.starsData);
            
            obj.currentNode = initNode;
            obj.currentTime = 0;
            obj.currentSettlementContext = SettlementContextEnum.FastShip;
            
            InitialObservation = obj.getObsInfoForTimeAndStarId(initTime, initStarId);
            
            % (optional) use notifyEnvUpdated to signal that the 
            % environment has been updated (e.g. to update visualization)
            notifyEnvUpdated(obj);
        end
        
        function obsInfo = getObsInfoForTimeAndStarId(obj, time, starId)                       
            [rVectKpc, vVectKpcMyr] = getStarPositionKpcMyr(starId, time, obj.starsData.starsData);
            star = obj.starsData.getStarForId(starId);
            
            obsInfo = {time, rVectKpc, vVectKpcMyr, [star.R; star.thetaF]};
%             obsInfo = {time, [star.R; star.thetaF]};
            
%             obsInfoBaseline = obj.getObservationInfo();
%             lb = [obsInfoBaseline.LowerLimit];
%             ub = [obsInfoBaseline.UpperLimit];
%             
%             for(i=1:length(obsInfo))
%                 subObsInfo = obsInfo{i};
%                 
%                 if(any(subObsInfo < lb(i)) || any(subObsInfo > ub(i)))
%                     a=1;
%                 end
%             end
            
%             obsInfo2 = obj.settlementGraph.getGalaxyStarObservationInfo();
%             obsInfo = horzcat(obsInfo1, num2cell(obsInfo2(:,4)'));
        end
    end

    methods(Static, Access=private)       
        function ObservationInfo = getInitObjInfo(settlementGraph)           
            %Time
            ObservationInfoTime = rlNumericSpec([1 1], 'LowerLimit',-3, 'UpperLimit',90);
            ObservationInfoTime.Name = 'Current Time [MYr]';
            ObservationInfoTime.Description = 'Current Time [Myr]';
            ObservationInfo = ObservationInfoTime;
            
            ObservationInfoCurrentPos = rlNumericSpec([3 1], 'LowerLimit',-32, 'UpperLimit',32);
            ObservationInfoCurrentPos.Name = 'Current Position [kpc]';
            ObservationInfoCurrentPos.Description = 'x, y, z [kpc]';
            ObservationInfo(end+1) = ObservationInfoCurrentPos;
            
            ObservationInfoCurrentVel = rlNumericSpec([3 1], 'LowerLimit',-1, 'UpperLimit',1);
            ObservationInfoCurrentVel.Name = 'Current Velocity [kpc/Myr]';
            ObservationInfoCurrentVel.Description = 'vx, vy, vz [kpc/Myr]';
            ObservationInfo(end+1) = ObservationInfoCurrentVel;
            
            ObservationInfoCurrentStarRThetaF = rlNumericSpec([2 1], 'LowerLimit',[2;-180], 'UpperLimit',[32;180]);
            ObservationInfoCurrentStarRThetaF.Name = 'Current Star R/Theta_F [kpc,deg]';
            ObservationInfoCurrentStarRThetaF.Description = 'R, Theta_F [kpc,deg]';
            ObservationInfo(end+1) = ObservationInfoCurrentStarRThetaF;
            
%             obsInfo = settlementGraph.getGalaxyStarObservationInfo();
%             ObservationInfoSettlement = rl.util.rlFiniteSetSpec.empty(1,0);
%             for(i=1:length(obsInfo(:,4)))
%                 id = obsInfo(i,1);
%                 fprintf('Creating observation for star %u\n', id);
%                 
%                 ObservationInfoSettlement(i) = rlFiniteSetSpec([0, 1]);
%                 ObservationInfoSettlement(i).Name = sprintf('Star %u Status', id);
%                 ObservationInfoSettlement(i).Description = '1 = settled star, 0 = not settled star';
%             end
%             ObservationInfo = [ObservationInfo, ObservationInfoSettlement];
        end
        
        function ActionInfo = getInitActionInfo(settlementGraph)
%             obsInfo = settlementGraph.getGalaxyStarObservationInfo();
%             ids = obsInfo(2:end,1);
            
            lb = [2, -180];
            ub = [32, 180];
            ActionInfo = rlNumericSpec([1 2], 'LowerLimit',lb, 'UpperLimit',ub);
            
%             ActionInfoStarSelection = rlFiniteSetSpec(ids);
%             ActionInfoStarSelection.Name = 'Choice of Stars to Settle';
%             ActionInfo = ActionInfoStarSelection;
%             
%             ActionInfoWaitDuration = rlNumericSpec([1 1], 'LowerLimit',3, 'UpperLimit',13);
%             ActionInfoWaitDuration.Name = 'Departure Wait Duration';
%             ActionInfoWaitDuration.Description = 'How long to wait before departing star?';
%             ActionInfo(end+1) = ActionInfoWaitDuration;
%             
%             ActionInfoFlightTime = rlNumericSpec([1 1], 'LowerLimit',0.1, 'UpperLimit',20);
%             ActionInfoFlightTime.Name = 'Flight Time Between Stars';
%             ActionInfoFlightTime.Description = 'How long is the duration of flight between stars?';
%             ActionInfo(end+1) = ActionInfoFlightTime;
%             
%             ActionInfoLambertTm = rlFiniteSetSpec([-1 1]);
%             ActionInfoLambertTm.Name = 'Lambert TM';
%             ActionInfo(end+1) = ActionInfoLambertTm;
        end
    end
end
