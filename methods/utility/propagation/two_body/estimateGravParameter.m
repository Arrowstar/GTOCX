function estMu = estimateGravParameter(rVects)   
	rMags = rssq(rVects,1);
    VCs = computeVc(rMags);
    MUs = (VCs.^2) .* rMags;
%     estMu = mean(MUs);    
    estMu = sum(MUs)/numel(MUs);
end