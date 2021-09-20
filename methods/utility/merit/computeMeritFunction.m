function [J,Er,Eth] = computeMeritFunction(ids, starData, actualDvs, maxDvs)
%Computes the merit function according to (26).  Only the middle factor,
%does not handle the delta-v ratio or the B time factor.
%   INPUTS:
%       ids - a list of star IDs that have been settled.
%       starData - the matrix of data from the stars.txt file
%       actualDvs - a Nx1 vector of the actual delta-v magnitudes executed
%           by every vehicle.  Each element is the total delta-v expended
%           for that vehicle.  Must correspond to the maxDvs array.
%       maxDvs - a Nx1 vector of the theoretical maximum delta-v magnitude 
%           executable by every vehicle.  Must correspond to the actualDvs array.
%   OUTPUTS:
%       J - the merit function value
%       Er - error function, R
%       Eth - error function, theta
    if(isempty(ids))
        J = 0;
        Er = Inf;
        Eth = Inf;
        return;
    end

    ids = unique(ids,'stable');

%     rowNum = ids+1;
%     rowData = starData(rowNum,:);

    [~,Locb] = ismember(ids,starData(:,1));
    rowData = starData(Locb,:);

    r_iArr = rowData(:,2);
    theta_iArr = deg2rad(rowData(:,6));

    N = length(r_iArr);
    
    Er = computeEr(r_iArr, N);
    Eth = computeEtheta(theta_iArr, N);
    
%     tEnd = datetime(2019,06,12,20,0,0);
%     tStart = datetime(2019,05,21,20,0,0);
%     tSub = datetime(); %take current time as submission time, close enough

%     tEnd = datenum(2019,06,12,20,0,0);
%     tStart = datenum(2019,05,21,20,0,0);
%     tSub = now();
%     B = 1 + ((tEnd - tSub)/(tEnd - tStart))^4;
    B = 1;
    
    deltaVRatio = sum(maxDvs)/sum(actualDvs);
    
    J = B * (N/(1 + (1E-4)*N*(Er+Eth)))*deltaVRatio;
end

function Er = computeEr(r_iArr, N)
    sum = 0;
    
    for(k=0:30)
        Rk = computeRk(k);
        sum = sum + (computeFr(Rk, r_iArr, N) / computeGr(Rk) - 1)^2;
    end
    
    Er = sum;
end

function Rk = computeRk(k)
    Rk = k + 2;
end

function fr = computeFr(r, r_iArr, N)
    sum = 0;
    sr = 1.0;
    
    for(i=1:N)
        r_i = r_iArr(i);
        sum = sum + computeF(r, r_i, sr);
    end
    
    fr = sum/N;
end

function gr = computeGr(r)
    Rmax = 32;
    Rmin = 2;

    gr = computeAlpha(r) * (2*r/(Rmax^2 - Rmin^2));
end

function alpha = computeAlpha(r)
    if(r == 2)
        alpha = 0.5833;
    elseif(r == 32)
        alpha = 0.4948;
    else
        alpha = 1;
    end
end

function Eth = computeEtheta(theta_iArr, N)
    sum = 0;
    
    for(k=0:32)
        thetaK = computeThetaK(k);
        sum = sum + (computeFtheta(thetaK, theta_iArr, N) / computeGtheta(thetaK) - 1)^2;
    end
    
    Eth = sum;
end

function thetaK = computeThetaK(k)
    thetaK = -pi + 2*pi*(k/32);
end

function ftheta = computeFtheta(theta,theta_iArr,N)
    sum = 0;
    stheta = 2*pi/32;
    
    for(i=1:N)
        theta_i = theta_iArr(i);
        sum = sum + computeF(theta, theta_i, stheta);
    end
    
    ftheta = sum/N;
end

function gtheta = computeGtheta(theta)
	gtheta = computeBeta(theta)/(2*pi);
end

function beta = computeBeta(theta)
    if(theta == -pi || theta == pi)
        beta = 0.5;
    else
        beta = 1;
    end
end

function f = computeF(x, mu,s)
    if(abs(x-mu)>=s)
        f = 0;
    else
        f = 1/s - abs(x-mu)/s^2;
    end
end