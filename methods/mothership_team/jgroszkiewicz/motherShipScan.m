% Set Scan Parameters
% Parameters are set for a coarse scan to test 
MinDeltaVMag = 50; % km/sec
MaxDeltaVMag = 500; % km/sec
StepDeltaVMag = 50; % km/sec
MinDeltaVTime = 0; % Myr
MaxDeltaVTime = 10; % Myr
StepDeltaVTime = 1; % Myr
MinDeltaVAzimuth = -180; % Degrees
MaxDeltaVAzimuth = 179.9999; % Degrees
StepDeltaVAzimuth = 30; % Degrees
MinDeltaVElevation = -30; % Degrees
MaxDeltaVElevation = 30; % Degrees
StepDeltaVElevation = 10; % Degrees
% Only evaluate stars whose radius from galactic center is within Rthr of
% the ship's current radius from galactic center
Rthr = 0.01; % kpc
% Constants
deg2rad = pi/180;
endTime = 90; % Myr
%
magArray = [MinDeltaVMag:StepDeltaVMag:MaxDeltaVMag]; nMag = length(magArray);
timeArray = [MinDeltaVTime:StepDeltaVTime:MaxDeltaVTime]; nTime = length(timeArray);
azArray = [MinDeltaVAzimuth:StepDeltaVAzimuth:MaxDeltaVAzimuth]; nAz = length(azArray);
elArray = [MinDeltaVElevation:StepDeltaVElevation:MaxDeltaVElevation]; nEl = length(elArray);
%
nScanPoints = nMag * nTime * nAz * nEl
shipId = zeros(1,nScanPoints);
shipDVxyz = zeros(3,nScanPoints);
shipDVtime = zeros(1,nScanPoints);
shipNEvalTimes = zeros(1,nScanPoints);
shipNStarsEval = zeros(1,nScanPoints);
shipValue = zeros(1,nScanPoints);
%
[solPos, solVel] = getStarPositionKpcMyr([0],timeArray',rawStarData);
%
nPointsScanned = 0;
timeScanStart = now;
for magIdx = 1:nMag
    mag = magArray(magIdx);
    for timeIdx = 1:nTime
        T = timeArray(timeIdx);
        shipPos0 = solPos(:,timeIdx);
        shipVel0 = solVel(:,timeIdx);
        for azIdx = 1:nAz
            az = azArray(azIdx)*deg2rad;
            sin_az = sin(az);
            cos_az = cos(az);
            for elIdx = 1:nEl
                el = elArray(elIdx)*deg2rad;
                sin_el = sin(el);
                cos_el = cos(el);
                %
                shipIdx = nPointsScanned+1;
                shipId(shipIdx) = shipIdx;
                shipDVtime(shipIdx) = T;
                shipDVxyz(:,shipIdx) = mag*[cos_az*cos_el; sin_az*cos_el; sin_el];
                shipVel0 = shipVel0 + shipDVxyz(:,shipIdx);
                %
                [evalTimes, evalStates, junk1, junk2, junk3] = propagateBody(shipDVtime(shipIdx),shipPos0,shipVel0,endTime-shipDVtime(shipIdx),[]);
                shipNEvalTimes(1,shipIdx) = size(evalTimes,1);
                %
                for evalIdx = 1:shipNEvalTimes(1,shipIdx)
                    T_eval = evalTimes(evalIdx,1);
                    Pos_eval = evalStates(evalIdx,1:3);
                    Vel_eval = evalStates(evalIdx,4:6);
                    R_eval = sqrt(Pos_eval*Pos_eval');
                end
                %
                nPointsScanned = nPointsScanned + 1;
                pctComplete = (nPointsScanned/nScanPoints)*100;
            end
        end
        disp(sprintf('%6.2f\n',pctComplete)); %#ok<DSPS>
    end
end
timeScanEnd = now;
scanDurationSec = (timeScanEnd - timeScanStart)*(86400);
disp(sprintf('Scan Duration (sec) = %15.2f',scanDurationSec)); %#ok<*DSPS>
