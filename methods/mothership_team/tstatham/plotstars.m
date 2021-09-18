% clear all
% close all
% clc

starData = readmatrix('data/stars.txt');

ids = starData(:,1);
count = 1;

for dmyr = 0:5:20
    tMyr = dmyr * ones(size(ids));
    [rVectKm, vVectKms] = getStarPositionKpcMyr(ids, tMyr, starData);
    

    subplot(2,2,count)
    plot(rVectKm(1,:), rVectKm(2,:), 'b.')
    %plot3(rVectKm(1,:), rVectKm(2,:), rVectKm(3,:), 'b.')

    title(sprintf('%i', dmyr)); 
    %title(fprintf('%d', count)); 
    count = count  + 1;
end


%thetaf = atan2d(rVectKm(2,:),rVectKm(1,:));
%err = abs(starData(ids+1,6) - thetaf');
%check = all(err < 1E-3);
%assert(check);