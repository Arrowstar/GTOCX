function [f, msCsvMatrix] = initShipTrajObjFcnFine(x, starIdList, starData)
    xx = [reshape(x',[2,numel(x)/2])',starIdList(:)];
    xx = reshape(xx',1,numel(xx));
    [~, ~, msDvs,spDvs, msCsvMatrix] = initShipTrajRunner(xx, true, starData);
    
    msDvMags = sum(msDvs,1);
    
    actualDvs = [msDvMags, spDvs(:)']';
%     maxDvs = [500*kmS2KpcMyr()*ones(1,size(allMs,3)), 300*kmS2KpcMyr()*ones([1,length(spDvs(:))])];
%     f = -1*computeMeritFunction(colonizedStarsIds, starData, actualDvs, maxDvs);
    f = sum(actualDvs); %#ok<UDIM>
end