classdef SettlementNode < matlab.mixin.SetGet & matlab.mixin.CustomDisplay
    %SettlementNode Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        star(1,1) Star
        settledOn(1,1) double %MYr
        
        maxNumAllowedDepartures(1,1) = 3
        departingSettlement(1,:) Settlement
        remainingSettlements(1,1) double = 3;
    end
    
    methods
        function obj = SettlementNode(star, settledOn, maxNumAllowedDepartures)
            if(nargin > 0)
                obj.star = star;
                obj.settledOn = settledOn;
                obj.maxNumAllowedDepartures = maxNumAllowedDepartures;
            end
        end
        
        function addDepartingSettlement(obj, settlement)
            obj.departingSettlement(end+1) = settlement;
            obj.remainingSettlements = obj.remainingSettlements - 1;
        end
    end
    
    methods(Access=protected)
        function displayScalarObject(obj)
            fprintf('Settlement of Star %u on %0.3f MYr\n', obj.star.id, obj.settledOn);
        end
    end
end