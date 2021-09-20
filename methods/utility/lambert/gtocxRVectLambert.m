function [vVect1, vVect2, trajOutsideBounds] = gtocxRVectLambert(rVect1, t1, rVect2, t2, tm)
%gtocxStarLambert Solves "Lambert's problem" between two position vectors in the Galaxy.
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

    dt = t2 - t1;

    [vVect1, ~] = gtocxRVectKeplerLambert(rVect1, t1, rVect2, t2, tm);
    
    fun = @(x) optObjFcn(x, t1, rVect1, dt, rVect2);
    x0 = vVect1;
    options = optimoptions(@fmincon, 'Display','off','Algorithm','interior-point','OptimalityTolerance',1E-12,'StepTolerance',1E-12);
	
    try
        [x, fval, exitflag,output] = fmincon(fun,x0,[],[],[],[],[],[],[],options);
        if(fval > 1E-4)
            warning('outside of star rend. limits: %0.9f\n%s', fval,output.message);
        end
    catch ME
        vVect1 = nan(size(rVect1));
        vVect2 = nan(size(rVect1));
        
        return
    end
    
    vVect1 = x(:);
    [~,y, ~,~,~] = propagateBody(t1,rVect1,vVect1,dt, [], []);
    
    radii = rssq(y(:,1:3),2);
    trajOutsideBounds = false;
    if(any(radii < 2 | radii > 32))
        trajOutsideBounds = true;
        warning('A trajectory was computed that has a position less than 2 kpc or greater than 32 kpc.');
    end
    
    vVect2 = y(end,4:6)';
end

function f = optObjFcn(x, t1, rVect1, dt, rVectStar2)
    vVect0 = x(1:3);
    vVect0 = vVect0(:);
    
    [t,y, ~,~,~] = propagateBody(t1,rVect1,vVect0,dt, [], []);
    
    rVectGuess = y(end,1:3)';
    
    f = norm(rVectStar2 - rVectGuess);
end