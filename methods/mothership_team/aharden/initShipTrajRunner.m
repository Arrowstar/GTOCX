function [allMs, colonizedStarsIds, msDvs,spDvs, msCsvMatrix] = initShipTrajRunner(x, useRealEoM, starData)
    numSettlePods = length(x)/(3*3);
    
    allMs = zeros(numSettlePods,3,3);
    allMs(:,:,1) = reshape(x(0*numSettlePods+1:3*numSettlePods),3,numSettlePods)';
    allMs(:,:,2) = reshape(x(3*numSettlePods+1:6*numSettlePods),3,numSettlePods)';
    allMs(:,:,3) = reshape(x(6*numSettlePods+1:9*numSettlePods),3,numSettlePods)';
    
    for(i=1:size(allMs,3))
        allMs(:,3,i) = starData(round(allMs(:,3,i)+1),1);
    end
    
    colonizedStarsIds = allMs(:,3,:);
    colonizedStarsIds = colonizedStarsIds(:);
    
    msDvs = zeros(numSettlePods,3);
    spDvs = zeros(numSettlePods,3);
    msDvVects = nan(3,numSettlePods,3);
    msDvTimes = nan(numSettlePods,3);
    spDvVects = nan(3,numSettlePods,3);
    spDvTimes = nan(numSettlePods,3);
    spStarIds = nan(numSettlePods,3);
    for(i=1:size(allMs,3))
        ms = allMs(:,:,i);
        
        starId0 = 0;
        t0 = ms(1,1);
        [rVect0t, vVect0t] = getStarPositionKpcMyr(starId0, t0, starData);
        
        for(j=1:size(ms,1))  
            starId1 = ms(j,3); 
            
            if(starId1 == starId0)
                msDvs(j,i) = NaN;
                spDvs(j,i) = NaN;
            else
                if(j == 1) %leaving from Earth
                    dt = ms(j,2); %only one burn segment here

                    [msDv, spDv, rVect2t, vVect2t, msDvVect, spDvVect] = computeBurnSegment(t0, dt, rVect0t, vVect0t, starId1, useRealEoM, starData);
                    
                    msDvs(j,i) = msDv;
                    spDvs(j,i) = spDv;
                    msDvVects(:,j,i) = msDvVect;
                    spDvVects(:,j,i) = spDvVect;
                    msDvTimes(j,i) = t0;
                    spDvTimes(j,i) = t0+dt;
                    spStarIds(j,i) = starId1;

                    segTotalDt = dt;
                else
                    dtCoast = ms(j,1);
                    dtBurn = ms(j,2);  
                    
                    [rVect1t, vVect1t] = computeCoastSegment(rVect0t, vVect0t, dtCoast, useRealEoM);
                    t1 = t0 + dtCoast;
                    [msDv, spDv, rVect2t, vVect2t, msDvVect, spDvVect] = computeBurnSegment(t1, dtBurn, rVect1t, vVect1t, starId1, useRealEoM, starData);
                    
                    msDvs(j,i) = msDv;
                    spDvs(j,i) = spDv;
                    msDvVects(:,j,i) = msDvVect;
                    spDvVects(:,j,i) = spDvVect;
                    msDvTimes(j,i) = t1;
                    spDvTimes(j,i) = t0+dtCoast+dtBurn;
                    spStarIds(j,i) = starId1;
                    
                    segTotalDt = dtCoast + dtBurn;
                end
                
                starId0 = starId1;
                t0 = t0+segTotalDt;
                rVect0t = rVect2t;
                vVect0t = vVect2t;
            end
        end
    end 
    
    msCsvMatrix = generateMotherShipCsvMatrix(numSettlePods, msDvTimes, spDvTimes, msDvVects, spDvVects, spStarIds);
    
    msDvs(isnan(msDvs)) = max(max(msDvs));
    spDvs(isnan(spDvs)) = max(max(spDvs));
end

function [msDv, spDv, rVect1t, vVect1t, msDvVect, spDvVect] = computeBurnSegment(t0, dt, rVect0t, vVect0t, starId1, useRealEoM, starData)
    [rVect1t, vVect1] = getStarPositionKpcMyr(starId1, t0+dt, starData);
    
    tm = 1;
    if(useRealEoM)
%         [~, vVect0tS, ~, vVect1tS] = gtocxStarLambert(starId0, t0, starId1, t0+dt, tm, starData);
        [vVect0tS, vVect1tS] = gtocxRVectLambert(rVect0t, t0, rVect1t, t0+dt, tm);
    else
%         [~, vVect0tS, ~, vVect1tS] = gtocxStarKeplerLambert(starId0, t0, starId1, t0+dt, tm, starData);
        [vVect0tS, vVect1tS] = gtocxRVectKeplerLambert(rVect0t, t0, rVect1t, t0+dt, tm);
    end
    dvShortVect = vVect0tS - vVect0t;
    dvShort = norm(dvShortVect);

    tm = -1;
    if(useRealEoM)
%         [~, vVect0tL, ~, vVect1tL] = gtocxStarLambert(starId0, t0, starId1, t0+dt, tm, starData);
        [vVect0tL, vVect1tL] = gtocxRVectLambert(rVect0t, t0, rVect1t, t0+dt, tm);
    else
%         [~, vVect0tL, ~, vVect1tL] = gtocxStarKeplerLambert(starId0, t0, starId1, t0+dt, tm, starData);
        [vVect0tL, vVect1tL] = gtocxRVectKeplerLambert(rVect0t, t0, rVect1t, t0+dt, tm);
    end
    dvLongVect = vVect0tL - vVect0t;
    dvLong = norm(dvLongVect);
    
    if(dvShort < dvLong)
        msDv = dvShort;
        spDvVect = vVect1 - vVect1tS;
        spDv = norm(spDvVect);
        vVect1t = vVect1tS;
        msDvVect = dvShortVect;
    else
        msDv = dvLong;
        spDvVect = vVect1 - vVect1tL;
        spDv = norm(spDvVect);
        vVect1t = vVect1tL;
        msDvVect = dvLongVect;
    end
end

function [rVect1, vVect1] = computeCoastSegment(rVect0t, vVect0t, dtCoast, useRealEoM)
    if(useRealEoM)
        [~,y, ~,~,~] = propagateBody(0, rVect0t, vVect0t, dtCoast, [], []);
        
        rVect1 = y(end,1:3);
        rVect1 = rVect1(:);
        vVect1 = y(end,4:6);
        vVect1 = vVect1(:);
    else
        gmu = estimateGravParameter(rVect0t);
        [rVect1,vVect1] = propagateTwoBodyMotion(rVect0t, vVect0t, dtCoast, gmu);
    end
end