classdef Star < matlab.mixin.SetGet
    %STAR Summary of this class goes here
    
    properties
        id
        R
        inc
        Om
        phi
        thF
        
        validColTarget = true;
        isColonized = false;
        isOptStar = false;
        colonizedBy = Star.empty(1,0);
        colonizedTime = NaN;
        timeSettlerShipToColonizeThisColonyDeparted = NaN;
        colDepartDvVect = NaN(3,1);
        colArriveDvVect = NaN(3,1);
        maxSettlerDv = 400;
        
        galaxy
    end
    
    properties(Dependent)
        totalColDv
        starData
    end
    
    methods
        function obj = Star(id, R, inc, Om, phi, thF, galaxy)
            obj.id = id;
            obj.R = R;
            obj.inc = inc;
            obj.Om = Om;
            obj.phi = phi;
            obj.thF = thF;
            obj.galaxy = galaxy;
            
            if(id == 0)
                obj.validColTarget = false;
            end
        end
        
        function starData = get.starData(obj)
            starData = obj.galaxy.starData;
        end
        
        function [rVectKpc, vVectKpcMyr] = getPositionVelocityAtTime(tMyr)
            [rVectKpc, vVectKpcMyr] = getStarPositionKpcMyr(obj.id, tMyr, obj.starData);
        end
        
        function colonizeStarsFromThisStar(obj, fH, isMonteCarlo)
            earliestColTime = obj.colonizedTime + 2;
            fprintf('Computing colony options from star %u...\n', obj.id);    
            fprintf('\tMinimum departure time: %0.3f Myr\n', earliestColTime);  

            [validStars] = obj.galaxy.getValidColonizableStarsForAStar(obj,earliestColTime);
            validStarIds = [validStars.id];
            
            numValidStars = length(validStars);
            fprintf('\tFound %u valid stars to compute trajectories to.\n', numValidStars); 
            
            gStarData = obj.starData;
            starId0 = obj.id;

            tt = tic;
            [posColStars, posColTDep, posColTArr, ~, posColDist, posColSpeed, posColDv1Vect, posColDv2Vect, ~, ~] = ...
                    Star.computeTrajectoriesToValidStars(starId0, earliestColTime, validStarIds, gStarData);
            tDur = toc(tt);
            
            fprintf('\tColonizing run took %0.3f seconds.\n', tDur);
            
            if(any(not(isnan(posColStars))))
                if(isMonteCarlo)
                    Inds = find(not(isnan(posColTArr)));
                    I = unique(Inds(randi(length(Inds),1,3)));
                else
                    [B,I] = mink(posColTArr,3);
                    I(isnan(B)) = [];
                end
                
                starsIdsToColonizeFromHere = posColStars(I)';
                starsToColFromHere = obj.galaxy.getStarsWithIDs(starsIdsToColonizeFromHere);
                starsDepartTimes = posColTDep(I);
                starsColTimes = posColTArr(I);
                starsColDepDv1Vect = posColDv1Vect(:,I);
                starsColDepDv2Vect = posColDv2Vect(:,I);
                starsColDists = posColDist(I);
                starsColSpeeds = posColSpeed(I);
                
                for(i=1:length(starsIdsToColonizeFromHere))
                    starsToColFromHere(i).setColonized(obj, starsColTimes(i), starsDepartTimes(i), starsColDepDv1Vect(:,i), starsColDepDv2Vect(:,i), starsColDists(i), starsColSpeeds(i), fH);
                end
                
                obj.displayCurrentJ();
                
                for(i=1:length(starsToColFromHere))
                    starsToColFromHere(i).colonizeStarsFromThisStar(fH);
                end
            else
                fprintf('No valid colony options found.  Moving on to next colony...\n');
                obj.displayCurrentJ();
            end
            
            drawnow;
        end 
        
        function colonized = colonizeFromThisStarToGivenStar(obj, star1)
            earliestColTime = obj.colonizedTime + 2;
            gStarData = obj.starData;
            starId0 = obj.id;
            starId1 = star1.id;
            
            [posColStars, posColTDep, posColTArr, posColDv, posColDist, posColSpeed, posColDv1Vect, posColDv2Vect, posColRVect1, posColRVect2] = ...
                    computeTrajectoriesToValidStars(starId0, earliestColTime, starId1, gStarData);
                
            if(all(not(isnan(posColStars))))
                colonized = true;
            else
                colonized = false;
            end
        end
                
        function J = computeCurrentJ(obj)
            colStars = obj.galaxy.getColonizedStars();
            colStarIds = [colStars.id];
            usedDv = [colStars.totalColDv];
            maxDv = [colStars.maxSettlerDv];
            [J,~,~] = computeMeritFunction(colStarIds, obj.starData, usedDv, maxDv);
        end
        
        function displayCurrentJ(obj)
            J = obj.computeCurrentJ();
            fprintf('Current J (no B): %0.3f \n', J);
            fprintf('Stars Colonized: %u \n', length(colStarIds));
        end
        
        function setColonized(obj, colByStar, colonizedTime, departTime, colDepartDvVect, colArriveDvVect, dist, speed, fH)
            obj.isColonized = true;
            obj.colonizedBy = colByStar;
            obj.colonizedTime = colonizedTime;
            obj.timeSettlerShipToColonizeThisColonyDeparted = departTime;
            obj.colDepartDvVect = colDepartDvVect;
            obj.colArriveDvVect = colArriveDvVect;
            obj.maxSettlerDv = 400;
            
            obj.writeColonizationDataToSubmissinoFile(fH);
            colPedStr = obj.computeColonizationPedigree();
            
%             [rVects,~] = getStarPositionKpcMyr([colByStar.id, obj.id]', [departTime, colonizedTime]', obj.starData);
%             dist = sqrt(sum(diff(rVects,1,2).^2));
%             speed = dist/(colonizedTime-departTime);
            
            fprintf('##################################################\n');
            fprintf('Colonizing star # %u...! %s\n', obj.id, colPedStr);
            fprintf('Colonized by: %u\n', colByStar.id);
            fprintf('Settler Ship Departure Time: %0.3f Myr\n', departTime);
            fprintf('Settler Ship Arrival Time: %0.3f Myr\n', colonizedTime);
            fprintf('Settler Ship Departure DV: %0.3f km/s\n', norm(colDepartDvVect));
            fprintf('Settler Ship Departure DV Vector: [%0.3f %0.3f %0.3f] km/s\n', colDepartDvVect(1), colDepartDvVect(2), colDepartDvVect(3));
            fprintf('Settler Ship Arrival DV: %0.3f km/s\n', norm(colArriveDvVect));
            fprintf('Settler Ship Arrival DV Vector: [%0.3f %0.3f %0.3f] km/s\n', colArriveDvVect(1), colArriveDvVect(2), colArriveDvVect(3));
            fprintf('Settler Ship Total DV: %0.3f km/s\n', norm(colDepartDvVect) + norm(colArriveDvVect));
            fprintf('Settler Ship Travel Distance: %0.3f kpc\n', dist);
            fprintf('Settler Ship Avg. Linear Travel Speed: %0.3f kpc/Myr\n', speed)
            fprintf('##################################################\n');
        end
        
        function writeColonizationDataToSubmissinoFile(obj, fH)
            if(obj.isColonized &&  not(isempty(fH)))
                fprintf(fH,'%u, %u, %u, %0.12f, %0.12f, %0.12f, %0.12f, %0.12f, %0.12f, %0.12f, %0.12f\n', ...
                            [obj.colonizedBy.id, ...
                            obj.id, ...
                            2, ...
                            obj.timeSettlerShipToColonizeThisColonyDeparted, ...
                            obj.colonizedTime, ...
                            obj.colDepartDvVect(:)', ...
                            obj.colArriveDvVect(:)']);
            end
        end
        
        function value = get.totalColDv(obj)
            value = norm(obj.colDepartDvVect) + norm(obj.colArriveDvVect);
        end
        
        function [colStr, colList] = computeColonizationPedigree(obj)
            if(obj.isColonized)
                colList = obj.id;
                
                colFrom = obj.colonizedBy;
                while(not(isempty(colFrom)))
                    colList(end+1) = colFrom.id; %#ok<AGROW>
                    
                    colFrom = colFrom.colonizedBy;
                end
                
                colStr = sprintf('(%s)',strjoin(strtrim(cellstr(num2str(colList'))),' < '));
            else
                colStr = '';
            end
        end
    end
    
    methods(Static)
        function [posColStars, posColTDep, posColTArr, posColDv, posColDist, posColSpeed, posColDv1Vect, posColDv2Vect, posColRVect1, posColRVect2] = ...
                    computeTrajectoriesToValidStars(starId0, earliestColTime, validStarIds, gStarData)
            numValidStars = length(validStarIds);
            
            posColStars = NaN(length(validStarIds),1);
            posColTDep = NaN(length(validStarIds),1);
            posColTArr = NaN(length(validStarIds),1);
            posColDv = NaN(length(validStarIds),1);
            posColDist = NaN(length(validStarIds),1);
            posColSpeed = NaN(length(validStarIds),1);
            
            posColDv1Vect = NaN(3,length(validStarIds));
            posColDv2Vect = NaN(3,length(validStarIds));
            posColRVect1 = NaN(3,length(validStarIds));
            posColRVect2 = NaN(3,length(validStarIds));
                                   
            parfor(i=1:numValidStars)
                starId1 = validStarIds(i);            
                
                [problem, x0] = Star.getOptProblemStruct(starId0, starId1, earliestColTime, gStarData);
                
                try
                    ms = MultiStart('Display','none');
                    tpoints = CustomStartPointSet(x0);
                    [bestX,~,bestExitFlag,~] = ms.run(problem,tpoints);

                    if(bestExitFlag == 1)
                        x = bestX;
                        
                        useHiFi = true;
                        if(isMonteCarlo)
                            useHiFi = false;
                        end

                        [c, ~, tempPosColDv, tempPosColDv1Vect, tempPosColDv2Vect, tempPosColRVect1, tempPosColRVect2] = Star.colonizeNonlcon(x, starId0, starId1, useHiFi, gStarData);
                        
                        if(all(c <= 1E-3))
                            posColStars(i) = starId1;
                            posColTDep(i) = x(1); 
                            posColTArr(i) = x(1) + x(2);

                            posColDv(i) = tempPosColDv;
                            posColDv1Vect(:,i) = tempPosColDv1Vect;
                            posColDv2Vect(:,i) = tempPosColDv2Vect;
                            
                            posColRVect1(:,i) = tempPosColRVect1;
                            posColRVect2(:,i) = tempPosColRVect2;
                            
                            [rVectAtT1_Star0] = getStarPositionKpcMyr(starId0, posColTArr(i), gStarData);
                            [rVectAtT1_Star1] = getStarPositionKpcMyr(starId1, posColTArr(i), gStarData);
                            
                            posColDist(i) = norm(rVectAtT1_Star1(1:2) - rVectAtT1_Star0(1:2));
                            posColSpeed(i) = posColDist(i)/(posColTArr(i) - posColTDep(i));
                        end
                    end
                catch ME
                    continue
                end
            end
        end
        
        function arrivalTime = colonizeObjFcn(x)
            arrivalTime = x(1) + x(2);
        end
        
        function [c, ceq, dvTot, dv1Vect, dv2Vect, rVect1, rVect2] = colonizeNonlcon(x, star0Id, star1Id, isHighFidelity, starData)
            c = [];
            ceq = [];
            
            [dvTot, dv1, dv2, dv1Vect, dv2Vect, rVect1, rVect2, trajOutsideBounds] = Star.dvForOneShip(x, star0Id, star1Id, isHighFidelity, starData);
            if(star0Id == 0)
                c(1) = dv1 - 199;
                c(2) = dv2 - 299;
            else
                if(isHighFidelity)
                    c(1) = dv1 - 174;
                    c(2) = dv2 - 174;
                else
                    c(1) = dv1 - 174;
                    c(2) = dv2 - 174;
                end
            end
            c(3) = double(trajOutsideBounds);
            
            c = c(:);
        end
        
        function [dvTot, dv1, dv2, dv1Vect, dv2Vect, rVect1, rVect2, trajOutsideBounds] = dvForOneShip(x, starId0, starId1, isHighFidelity, starData)
            t0 = x(1);
            dt = x(2);
            t1 = t0+dt;

            [dvTot1, dv11, dv12, dv11Vect, dv12Vect, rVect11, rVect12, trajOutsideBounds1] = Star.computeLambertArcDv(starId0, t0, starId1, t1, 1, isHighFidelity, starData);
            [dvTot2, dv21, dv22, dv21Vect, dv22Vect, rVect21, rVect21, trajOutsideBounds2] = Star.computeLambertArcDv(starId0, t0, starId1, t1, -1, isHighFidelity, starData);

            if(dvTot1 < dvTot2)
                dvTot = dvTot1; 
                dv1 = dv11;
                dv2 = dv12;
                dv1Vect = dv11Vect;
                dv2Vect = dv12Vect;
                rVect1 = rVect11;
                rVect2 = rVect12;
                trajOutsideBounds = trajOutsideBounds1;
            else
                dvTot = dvTot2;
                dv1 = dv21;
                dv2 = dv22;
                dv1Vect = dv21Vect;
                dv2Vect = dv22Vect;
                rVect1 = rVect21;
                rVect2 = rVect22;
                trajOutsideBounds = trajOutsideBounds2;
            end
        end

        function [dvTot, dv1, dv2, dv1Vect, dv2Vect, rVect0, rVect1, trajOutsideBounds] = computeLambertArcDv(starId0, t0, starId1, t1, tm, isHighFidelity, starData)
            [rVect0, vVect0] = getStarPositionKpcMyr(starId0, t0, starData);
            [rVect1, vVect1] = getStarPositionKpcMyr(starId1, t1, starData);

            if(isHighFidelity)
                [~, vVect0t, ~, vVect1t, trajOutsideBounds] = gtocxStarLambert(starId0, t0, starId1, t1, tm, starData);
            else
                [~, vVect0t, ~, vVect1t] = gtocxStarKeplerLambert(starId0, t0, starId1, t1, tm, starData);
                trajOutsideBounds = false;
            end

            dv1Vect = (vVect0t-vVect0)*(1/kmS2KpcMyr);
            dv1 = norm(dv1Vect);
            
            dv2Vect = (vVect1-vVect1t)*(1/kmS2KpcMyr);
            dv2 = norm(dv2Vect);
            
            dvTot = (dv1 + dv2);
        end
        
        function [problem, x0] = getOptProblemStruct(starId0, starId1, earliestColTime, gStarData)
            fun = @(x) Star.colonizeObjFcn(x);
            options = optimoptions(@fmincon, 'Display','off', 'OptimalityTolerance',1E-2, 'StepTolerance',1E-2, 'ConstraintTolerance',1E-3);
            
            A = [1 1];
            b = 90;

            lb = [earliestColTime, 1];
            ub = [90, 90];
            nonlcon = @(x) Star.colonizeNonlcon(x, starId0, starId1, false, gStarData);

%                 xMid = (ub-lb).*rand(3,2) + lb;
%                 x0 = [lb; xMid; ub];
            x0 = [lb; ub];

            problem = createOptimProblem('fmincon','objective',fun, 'x0',lb, 'Aineq',A, 'bineq',b, 'lb',lb, 'ub',ub, 'nonlcon',nonlcon, 'options',options);
        end
    end
end

