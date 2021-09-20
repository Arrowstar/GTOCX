function [rVect1, vVect1, rVect2, vVect2, trajOutsideBounds] = gtocxStarLambert(starId1, tStar1, starId2, tStar2, tm, starData)
%gtocxStarLambert Solves "Lambert's problem" between stars in the Galaxy.
%   INPUTS:
%       starId1 - ID number of star 1, the departure star
%       tStar1 - time of departure from Star 1 in Myr
%       starId2 - ID number of star 2, the arrival star
%       tStar2 - time of arrival at Star 2 in Myr. Must be greater than
%               tStar1
%       tm - Lmabert like long way or short way.  Use +1 for short way and
%               -1 for long way.  See lambert.m.
%       starData - The matrix of data from the stars.txt file.
%   OUTPUTS:
%       rVect1 - The position vector of the vehicle at tStar1 (co-located
%           with Star 1).
%       vVect1 - The velocity vector of the vehicle at tStar1 such that it
%           will intersect with Star 2 in (tStar2 - tStar1) Myr.
%       rVect2 -  The position vector of the vehicle at tStar2 (co-located
%           with Star 2).     
%       vVect2 - The velocity vector of the vehicle at tStar2.

    dt = tStar2 - tStar1;

    [R1, V1, R2, ~] = gtocxStarKeplerLambert(starId1, tStar1, starId2, tStar2, tm, starData);
    
    rVect0 = R1;
    vVect0 = V1;
    
    fun = @(x) optObjFcn(x, tStar1, rVect0, dt, R2);
    x0 = vVect0;
    options = optimoptions(@fmincon, 'Display','off','Algorithm','interior-point','OptimalityTolerance',1E-12,'StepTolerance',1E-12);
	
    try
        [x, fval, exitflag, output] = fmincon(fun,x0,[],[],[],[],[],[],[],options);
        if(fval > 1E-4)
            warning('outside of star rend. limits: %0.9f\n%s', fval,output.message);
        end
    catch ME
        rVect1 = R1;
        rVect2 = R2;
        vVect1 = nan(size(rVect1));
        vVect2 = nan(size(rVect1));
        
        return
    end
    
    vVect0 = x(:);
    [~,y, ~,~,~] = propagateBody(tStar1,rVect0,vVect0,dt, [], []);
    
    radii = rssq(y(:,1:3),2);
    trajOutsideBounds = false;
    if(any(radii < 2 | radii > 32))
        trajOutsideBounds = true;
%         warning('A trajectory was computed that has a position less than 2 kpc or greater than 32 kpc.');
    end
    
    rVect1 = y(1,1:3)';
    rVect2 = y(end,1:3)';
    
    vVect1 = x(:);
    vVect2 = y(end,4:6)';
end

function f = optObjFcn(x, tStar1, rVect0, dt, rVectStar2)
    vVect0 = x(1:3);
    vVect0 = vVect0(:);
    
    [~,y, ~,~,~] = propagateBody(tStar1,rVect0,vVect0,dt, [], []);
    
    rVectGuess = y(end,1:3)';
    
    f = norm(rVectStar2 - rVectGuess);
end