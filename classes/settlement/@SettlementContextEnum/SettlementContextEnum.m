classdef SettlementContextEnum < matlab.mixin.SetGet
    %SettlementContextEnum Summary of this class goes here
    %   Detailed explanation goes here
    
    enumeration
        Mothership('Mothership', 200*kmS2KpcMyr(), 500*kmS2KpcMyr(), 3)
        SettlePod('Settlement Pod', 300*kmS2KpcMyr(), 300*kmS2KpcMyr(), 1)
        FastShip('Fast Ship', 1500*kmS2KpcMyr(), 1500*kmS2KpcMyr(), 2)
        SettlerShip('Settler Ship', 175*kmS2KpcMyr(), 400*kmS2KpcMyr(), 3);
    end
    
    properties
        name(1,:) char
        maxIndivDv(1,1) double
        maxTotalDv(1,1) double
        maxNumImpulses(1,1) double
    end
    
    methods
        function obj = SettlementContextEnum(name, maxIndivDv, maxTotalDv, maxNumImpulses)
            obj.name = name;
            obj.maxIndivDv = maxIndivDv;
            obj.maxTotalDv = maxTotalDv;
            obj.maxNumImpulses = maxNumImpulses;
        end
    end
end

