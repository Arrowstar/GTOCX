function [c,ceq] = initShipTrajNonlconCoarse(x, starData)
    [~,colonizedStarsIds,msDvs,spDvs] = initShipTrajRunner(x, false, starData);
        
    msDvSums = sum(msDvs,1);    
    msDvsVect = msDvs(:);
    numNonZeroDvs = sum(msDvs > 0.01*kmS2KpcMyr(),1);
    
    c = [];
    c = msDvsVect - 200*kmS2KpcMyr();
    c = vertcat(c, msDvSums(:) - 500*kmS2KpcMyr());
    c = vertcat(c, spDvs(:) - 300*kmS2KpcMyr());
    c = vertcat(c, 1000*(length(colonizedStarsIds) - length(unique(colonizedStarsIds))));
    
    if(max(x)==2)
        c  = 100000.*ones(size(c));
    end

    
    ceq = [];
end