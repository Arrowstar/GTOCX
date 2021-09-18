function [f, msCsvMatrix] = initShipTrajObjFcnCoarse(x, starData)
    [allMs, colonizedStarsIds, msDvs,spDvs, msCsvMatrix] = initShipTrajRunner(x, false, starData);
    
    msDvMags = sum(msDvs,1);
    
    actualDvs = [msDvMags, spDvs(:)']';
    maxDvs = [500*kmS2KpcMyr()*ones(1,size(allMs,3)), 300*kmS2KpcMyr()*ones([1,length(spDvs(:))])];
    f = -1*computeMeritFunction(colonizedStarsIds, starData, 1, 1);
%     f = sum(actualDvs); %#ok<UDIM>

    if(isnan(f) || not(isfinite(f)))
        error('F func not a number or not finite.');
    end
end