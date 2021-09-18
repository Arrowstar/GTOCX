clc; clear variables; clear format long g; close all; fclose('all');

% addpath(genpath('data'));
% addpath(genpath('methods'));
% addpath(genpath('solutions'));
% addpath(genpath('studies'));
addpath(genpath('.'));

starData = dlmread('data/stars.txt',',',1,0);
optStarIds = unique(csvread('data/solN2703.csv'));
% optStarIds = [];

galaxy = Galaxy(starData, optStarIds);

% initColIDs = [13412, 10450, 96537, 10525, 64521];
% colTimes = [17.704083747983, 44.749759475423, 35.987374282272, 34.594127482005, 35.213914928679]';

% initColIDs = [34592, 9489, ...
%              13412, 86403, ...
%              96537, 20154, ...
%              60505, 72094, 50379];
% colTimes = [25.743043761103, 32.199187017494, ...
%             14.850706490168, 30.270169178521, ...
%             14.330218862515, 31.203810479918, ...
%             27.635938520633, 53.135994083504, 65.336668990620]';

% initColIDs = [34592, 16275, ...
%              13412, 86403, ...
%              28620, 41971, 62024, ...
%              60505, 72094, 50379];
% colTimes = [25.743043761103, 48.762465237146, ...
%             14.850706490168, 30.270169178521, ...
%             15.9582265625, 30.1582, 53.4582,...
%             27.635938520633, 53.135994083504, 65.336668990620]';

% initColIDs = [35916,79923,952, ...
%               88169,9423,39333, ...
%               5169,21241,48243, ...
%               34592,16275];
% colTimes = [10.410455654338126, 20.410455654338126, 30.410455654338126,...
%             11.425829636055333, 21.425829636055333, 31.425829636055333,...
%             16.758661252981742, 26.758661252981742, 36.75866125298174,...
%             31.148333895925543, 26.160459029443874];

initColIDs = [11125;27230;43144;59172;27686];
colTimes = 15 * ones(size(initColIDs));

isMonteCarlo = true;

stars = galaxy.getStarsWithIDs(initColIDs);
[stars(:).isColonized] = deal(true);
[stars(:).colonizedBy] = deal(galaxy.getStarsWithIDs(0));

for(i=1:length(colTimes))
    stars(i).colonizedTime = colTimes(i);
end

[stars(:).colDepartDvVect] = deal([500;0;0]);
[stars(:).colArriveDvVect] = deal([500;0;0]);
[stars(:).maxSettlerDv] = deal(500);

[fH,filename] = galaxy.openSettlerShipOutputFile();
[~,name,~] = fileparts(filename);
diary(sprintf('%s.log',name));
fprintf('Starting colonization.  Writing results to file: %s\n', filename);
t = tic;
for(i=1:length(initColIDs))
    stars(i).colonizeStarsFromThisStar(fH, isMonteCarlo);
end
tDur = toc(t);
fclose(fH);
diary off;

save(sprintf('%s.mat',name));
fprintf('Total Run Time = %0.3f s\n', tDur);
