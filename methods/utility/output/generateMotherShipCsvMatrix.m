function msCsvMatrix = generateMotherShipCsvMatrix(numSPods, msDvTimes, spDvTimes, msDvVects, spDvVects, spStarIds)
% generateMotherShipCsvMatrix A function that generates the output CSV data
% for the mothership portion of the solution file.
%   INPUTS:
%       numSPods - number of settle pods per mothership
%       msDvTimes - NxM matrix of times of the mothership delta-vs, where N
%                   is the number of delta-vs and M is the number of
%                   motherships [Myr]
%       spDvTimes - NxM matrix of the times of the settlement pods
%                   delta-vs, where N is the Nth pod from the Mth
%                   mothership [Myr]
%       msDvVects - 3xNxM matrix of the delta-v vectors of the mothership
%                   burns, where N is the Nth burn of the Mth moterhship
%                   [kpc/Myr]
%       spDvVects - 3xNxM matrix of the delta-v vectors of the settle pod
%                   burns, where is the Nth settle pod from the Mth moterhship
%                   [kpc/Myr]
%       spStarIds - an NxM matroc of the IDs of stars settled by the settle
%                   pods, where N is the Nth settle pod from the Mth
%                   mothership [-]
%   OUTPUTS:
%       msCsvMatrix - a Nx1 cell array, with each cell containing a row
%                     vector that corresponds to a line in the solution
%                     file.  Use writecell() to write to a file.  
%                     EG: writecell(msCsvMatrix,'msOut20190521_3MS_2SP.txt');

    msCsvMatrix = {};
    for(i=1:size(msDvTimes,2))
        dvTimes = msDvTimes(:,i);
        dvVects = msDvVects(:,:,i);
        msCsvMatrix{end+1} = generateMsLine(i, numSPods, dvTimes, dvVects); %#ok<AGROW>
        
        for(j=1:numSPods)
            starId = spStarIds(j,i);
            dvTime = spDvTimes(j,i);
            dvVect = spDvVects(:,j,i);
            
            msCsvMatrix{end+1} = generateSpLine(i, j, starId, dvTime, dvVect); %#ok<AGROW>
        end
    end
    
    msCsvMatrix = msCsvMatrix(:);
end

function msRow = generateMsLine(msId, numSPods, dvTimes, dvVects)
    msRow = nan(1,4);
    
    M = numel(dvTimes);
    
    msRow(1) = -1*msId;
    msRow(2) = 0;
    msRow(3) = numSPods;
    msRow(4) = numel(dvTimes);
    msRow(5:5+M-1) = dvTimes(:)';
    
    for(i=1:M)
        msRow(end+1:end+3) = dvVects(:,i)' * (1/kmS2KpcMyr());
    end    
end

function spRow = generateSpLine(msId, spId, starId, dvTime, dvVect)
    spRow = nan(1,7);
    
    spRow(1) = -1*msId;
    spRow(2) = spId;
    spRow(3) = starId;
    spRow(4) = dvTime;
    spRow(5:7) = dvVect(:)' * (1/kmS2KpcMyr());
end
