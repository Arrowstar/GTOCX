function E = keplerEq(M,e,eps)
% Function solves Kepler's equation M = E-e*sin(E)
% Input - Mean anomaly M [rad] , Eccentricity e and Epsilon 
% Output  eccentric anomaly E [rad]. 
   	En  = M;
	Ens = En - (En-e*sin(En)- M)/(1 - e*cos(En));
    cnt = 0;
	while ( abs(Ens-En) > eps && cnt <= 500)
		En = Ens;
		Ens = En - (En - e*sin(En) - M)/(1 - e*cos(En));
        cnt = cnt + 1;
	end
	E = Ens;