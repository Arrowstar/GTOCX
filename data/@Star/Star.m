classdef Star < matlab.mixin.SetGet
    %Star Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        id(1,1) double
        R(1,1) double
        i(1,1) double
        Omega(1,1) double
        phi(1,1) double
        thetaF(1,1) double
    end
    
    methods
        function obj = Star(id, R, i, Omega, phi, thetaF)
            if(nargin > 0)
                num = length(id);
                obj(num) = obj;
                
                if(num > 1)
                    obj(1) = Star();
                end

                for(j=1:num)
                    obj(j).id = id(j);
                    obj(j).R = R(j);
                    obj(j).i = i(j);
                    obj(j).Omega = Omega(j);
                    obj(j).phi = phi(j);
                    obj(j).thetaF = thetaF(j);
                end
            end
        end
    end
end