function [full_results, submission_results] = fastship_maneuver(starIds, t0Vec, dmyrVec, tm, starData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
% starIds = target stars for the lambert solver. Should match stars.txt file
% t0Vec = vector of potential launch times. Can be a single number or a range
% dymr = vector of potential travel times. Can be a single number or a range
% tm = 1 or -1, short or long lambert
% starData = the data file submitted for the program
%
%OUTPUTS:
% submission results = results in the format for submission according to
% GTOC
% full_results = submission results + more interesting info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

starId0 = 0;
full_results = [];
submission_results = [];

parfor i = 1:1:length(starIds) 
    starId1 = starIds(i);
    for t0 = t0Vec
        
        [~, vVect0t] = getStarPositionKpcMyr(starId0, t0, starData);
        
        % Years to go through
        for dmyr = dmyrVec
        
            %Get destination star velocity
            [~, vVect1] = getStarPositionKpcMyr(starId1, t0 + dmyr, starData); %kpc
            
            %Use short method
            [~, vVect0tS, ~, vVect1tS] = gtocxStarLambert(starId0, t0, starId1, t0+dmyr , tm, starData); 
            
            %Short DeltaV to get to the destination
            dvi = norm(vVect0tS - vVect0t)*(30856775814671900)/1000000/31557600;
            v1xyz =  (vVect0tS - vVect0t)*(30856775814671900)/1000000/31557600;
            
            %DeltaV required to match the stars velocity
            dvf = norm(vVect1 - vVect1tS)*(30856775814671900)/1000000/31557600;
            v2xyz =  (vVect1 - vVect1tS)*(30856775814671900)/1000000/31557600;

            %Calculate the total delta V 
            delVtotal = (dvi + dvf); %km/s
        
            %%r = deg2rad(rowData(:,2));
            %%theta = deg2rad(rowData(:,6));
        
            %Create a table output
            full_results = [full_results;[t0, dmyr, starId1, v1xyz(1), v1xyz(2), v1xyz(3), v2xyz(1), v2xyz(2), v2xyz(3), delVtotal]];
            submission_results = [submission_results; [starId1, t0, t0 + dmyr, v1xyz', v2xyz']];
            %save(filename,'position_table','-append');
        
        end
    end
end
