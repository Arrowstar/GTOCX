function Vc = computeVc(r)
%Computes the circular speed of the star given its radius.
%   INPUTS:
%       r -     The radius of the star in kpc
%   OUTPUTS:
%       Vc -    The circular velocity, Vc, in kpc/Myr according to (1).  
%               Note that (1) actually produces velocities in km/s so this 
%               is converted to the correct units.  
    k0 = 0.00287729;
    k1 = 0.0023821;
    k2 = -0.0010625;
    k3 = 0.000198502;
    k4 = -1.88428e-05;
    k5 = 9.70521e-07;
    k6 = -2.70559e-08;
    k7 = 3.7516e-10;
    k8 = -1.94316e-12;
    
    Vc = (1./(k8.*r.^8 + ...
               k7.*r.^7 + ...
               k6.*r.^6 + ...
               k5.*r.^5 + ...
               k4.*r.^4 + ...
               k3.*r.^3 + ...
               k2.*r.^2 + ...
               k1.*r.^1 + ...
               k0)) * kmS2KpcMyr();
end