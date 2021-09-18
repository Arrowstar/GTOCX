function [c,ceq] = initShipTrajNonlconFine(x, starIdList, starData)
    xx = [reshape(x',[2,numel(x)/2])',starIdList(:)];
    xx = reshape(xx',1,numel(xx));
    
    [~,~,msDvs,spDvs] = initShipTrajRunner(xx, true, starData);
        
    msDvSums = sum(msDvs,1);    
    msDvsVect = msDvs(:);
    numNonZeroDvs = sum(msDvs > 0.01*kmS2KpcMyr(),1);
    
    c = msDvsVect - 200*kmS2KpcMyr();
    c = vertcat(c, msDvSums(:) - 500*kmS2KpcMyr());
    c = vertcat(c, spDvs(:) - 300*kmS2KpcMyr());
%     c = vertcat(c, numNonZeroDvs(:) - 3);
    
    ceq = [];
end