function [rVect1, vVect1, rVect2, vVect2] = gtocxStarKeplerLambert(starId1, tStar1, starId2, tStar2, tm, starData)
%gtocxStarLambert Solves the actual Keplerian Lambert's problem between 
%                 stars in the Galaxy.  Motion provided from this function
%                 is KEPLERIAN, not according to the problem's equations of
%                 motion.  Use this function to get a very fast rough 
%                 estimate of the trajectory without needing to use the
%                 real EOMs, ode45, and fmincon.
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
    persistent mlVer
    if(isempty(mlVer))
        mlVer = version('-release');
    end
    
	if(tm > 0)
        tm = 1;
    else
        tm = -1;
	end

    dt = tStar2 - tStar1;

    ids = [starId1; starId2];
    tMyr = [tStar1; tStar2];
    [rVects, ~] = getStarPositionKpcMyr(ids, tMyr, starData);

	estMu = estimateGravParameter(rVects);
    
    R1 = rVects(:,1)';
    R2 = rVects(:,2)';
    
    if(ispc())
        if(strcmpi(mlVer,'2017b'))
            [vVect1, vVect2] = lambert_mex_R2017b(R1,R2,tm*dt,0,estMu);
        elseif(strcmpi(mlVer,'2019a'))
            [vVect1, vVect2] = lambert_mex_R2019a(R1,R2,tm*dt,0,estMu);
        else
            [vVect1, vVect2] = lambert(R1,R2,tm*dt,0,estMu);
        end
    else
        [vVect1, vVect2] = lambert(R1,R2,tm*dt,0,estMu);
    end
    
    rVect1 = R1(:);
    rVect2 = R2(:);
    vVect1 = vVect1(:);
    vVect2 = vVect2(:);
end