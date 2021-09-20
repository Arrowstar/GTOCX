classdef StarsData < matlab.mixin.SetGet
    %StarsData Summary of this class goes here
    %   Detailed explanation goes here
    
    properties       
        stars(:,1) Star
        starsData(:,6) double
    end
    
    methods
        function obj = StarsData(starsDataLoc)
            starsData = dlmread(starsDataLoc, ',', 1,0);
            obj.stars = Star(starsData(:,1), starsData(:,2), starsData(:,3), starsData(:,4), starsData(:,5), starsData(:,6));
            
            obj.starsData = starsData;
        end
        
        function star = getStarForId(obj, id)
            star = obj.stars([obj.stars.id] == id);
        end
        
        function starIds = getStarClosestTo(obj, R, thetaF)
            allR = obj.starsData(:,2);
            rNormalized = (allR - min(allR))/(max(allR) - min(allR));
            normR = (R - min(allR))/(max(allR) - min(allR));
            
            allT = obj.starsData(:,6);
            tNormalized = (allT - min(allT))/(max(allT) - min(allT));
            normTf = (thetaF - min(allT))/(max(allT) - min(allT));
            
            allVects = [rNormalized(:)'; tNormalized(:)'];
            thisVect = [normR;normTf];
            
            distVect = allVects - thisVect;
            distances = sqrt(sum(distVect.^2, 1));
            [~,I] = sort(distances);
            starIds = obj.starsData(I,1);
        end
        
        function obsInfo = getTemplateObsMatrix(obj)
            obsInfo = obj.starsData(:,[1,2,6]);
            obsInfo = [obsInfo, zeros(height(obsInfo),1)];
        end
    end
end

