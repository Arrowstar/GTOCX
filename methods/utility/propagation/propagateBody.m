function [t,y, te,ye,starIds] = propagateBody(t0,rVect0,vVect0,dt, step, starData)
%propagateBody Propagates a body according to the equations of motion in
%   (3)-(5).
%   INPUTS:
%       t0 - The time, in Myr, that the initial state is provided at.
%       Relative to Year 0.
%       rVect0 - The initial 3x1 position vector, in kpc
%       vVect0 - The initial 3x1 velocity vector, in kpc/Myr
%       dt - The duration of time to propagate for, in Myr
%       step - The step size to return the propagation for, in Myr.  If
%       empty, the solver just returns the native ode45 output.
%       starData - (OPTIONAL) The matrix of data from the stars.txt file.
%       If provided, this function will attempt to find star rendezvous
%       times and states, if they exist.  If this input is not needed, just
%       input an empty matrix [].  If empty, the propagation will not find
%       star rendezvous.
%   OUTPUTS:
%       t - Time of each state output from the  propagation, in Myr.
%       y - State of each output time from the propagation, in kpc and
%       kpc/Myr.  Each row is a state.  The first three elements are
%       position, the last three are velocity.
%       te - If starData is not empty, the times of each rendezvous event.
%       May be empty if no rendezvous found or starData is empty.
%       ye - The state of the vehicle at each rendezvous event, given as a
%       Nx6 matrix.  May be empty as above.
%       starIds - IDs of the stars found in rendezvous events.  May also be
%       empty as above.

    persistent odefun optionsBase
    if(isempty(odefun))
        odefun = @(t,y) gtocxOdefun(t,y);
        optionsBase = odeset('RelTol',1E-9, 'AbsTol',1E-9, 'InitialStep',dt, 'MaxStep',Inf, 'Refine',1);
    end

    if(dt == 0)
        t = t0;
        y = [rVect0', vVect0'];
        
        te = [];
        ye = [];
        starIds = [];
        
        return;
    end
    
    options = optionsBase;
    if(not(isempty(starData)))
        ids = starData(:,1);
        options = odeset(options, 'Events', @(t,y) starRendezvousEvts(t,y, starData,ids));
    end    
    
    if(isempty(step))
        tspan = [t0 t0+dt];
    else
        tspan = [t0:step:t0+dt];
    end
    y0 = [rVect0;vVect0];
        
    if(not(isempty(starData)))
        %call the odefun WITH the events finder code
        [t,y, te,ye,ie] = ode45(odefun,tspan,y0,options);
        starIds = ids(ie);
    else
        %call the odefun WITHOUT the events finder code
        [t,y] = ode45(odefun,tspan,y0,options);
        te = [];
        ye = [];
        starIds = [];
    end
end

function [value,isterminal,direction] = starRendezvousEvts(t,yvect, starData,ids)
    rVect = yvect(1:3);
    rVect = rVect(:);
    
    [rVectStars, ~] = getStarPositionKpcMyr(ids, t, starData);
    dist=rssq(rVect - rVectStars,1) - 1E-4;
    
    value = dist;
    isterminal = zeros(size(value));
    direction = zeros(size(value));
    
%     fprintf('%0.6f - %0.6f\n',t,min(dist));
end
