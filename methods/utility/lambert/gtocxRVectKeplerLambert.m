function [vVect1, vVect2] = gtocxRVectKeplerLambert(rVect1, t1, rVect2, t2, tm)
%gtocxStarLambert Solves the actual Keplerian Lambert's problem between 
%                 position vectors in the Galaxy.  Motion provided from this function
%                 is KEPLERIAN, not according to the problem's equations of
%                 motion.  Use this function to get a very fast rough 
%                 estimate of the trajectory without needing to use the
%                 real EOMs, ode45, and fmincon.
%   INPUTS:
%       rVect1 - 3x1 position vector of initial point, kpc
%       t1 - the time at position 1, Myr
%       rVect2 - 3x1 position vector of final point, kpc
%       t2 - the time at position 2, Myr. Must be greater than t1
%       tm - Lmabert like long way or short way.  Use +1 for short way and
%               -1 for long way.  See lambert.m.
%   OUTPUTS:
%       vVect1 - The velocity vector of the vehicle at t1 such that it
%           will intersect with rVect2 in (t2 - t1) Myr.   
%       vVect2 - The velocity vector of the vehicle at t2.
    persistent mlVer
    if(isempty(mlVer))
        mlVer = version('-release');
    end
    
	if(tm > 0)
        tm = 1;
    else
        tm = -1;
	end

    dt = t2 - t1;

    rVects = [rVect1, rVect2];
    
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
    
    vVect1 = vVect1(:);
    vVect2 = vVect2(:);
end