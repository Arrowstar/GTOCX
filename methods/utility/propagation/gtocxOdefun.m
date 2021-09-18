function [ydot] = gtocxOdefun(t,yvect)
%The function used by ODE45 to propagate a body.  See equations (3)-(5).
    %INPUTS: t - time in Myr
    %        yvect - 6-vector of position in velocity in kpc and kpc/Myr
    %OUTPUTS:
    %       ydot - the derivatives of the equations of motion in equations
    %       (3)-(5).

    x = yvect(1);
    y = yvect(2);
    z = yvect(3);
    vx = yvect(4);
    vy = yvect(5);
    vz = yvect(6);
    
    r = sqrt(x^2 + y^2 + z^2);
    vc = computeVc(r);
    f = vc^2/r; 
    
    ydot = zeros(6,1);
    
    ydot(1) = vx;
    ydot(2) = vy;
    ydot(3) = vz;
    
    ydot(4) = (-x/r)*f; %(3)
    ydot(5) = (-y/r)*f; %(4)
    ydot(6) = (-z/r)*f; %(5)
    
    ydot = ydot(:);
end