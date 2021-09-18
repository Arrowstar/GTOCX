function [rVect1,vVect1] = propagateTwoBodyMotion(rVect0, vVect0, dt, gmu)
   
%     [sma, ecc, inc, raan, arg, tru0] = getKeplerFromState(rVect0,vVect0,gmu);
%     mean0 = computeMeanFromTrueAnom(tru0, ecc);
%     meanMotion = computeMeanMotion(sma, gmu);
%     
%     mean1 = mean0 + meanMotion.*dt;
%     
%     tru1 = computeTrueAnomFromMean(mean1, ecc);
%     [rVect1,vVect1]=vect_getStatefromKepler(sma, ecc, inc, raan, arg, tru1, gmu);

    [rVect1,vVect1] = keplerUniversal(rVect0, vVect0, dt, gmu);
end