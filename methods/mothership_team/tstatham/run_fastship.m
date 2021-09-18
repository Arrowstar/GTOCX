clear all
close all
clc

%START TIMER
tic;
 
starData = readmatrix('data/stars.txt');  %Read in Stardata 
ids = starData(:,1);                      %starids
starId0 = 0;                              %home star, sol
%targetIDS = [50162, 54172, 10239,55720, 56277, 79921, 48198, 84052, 37153, 83446, 65846, 55720, 77256, 13634, 5397, 87651, 51919];                %Target Stars: single number or vector

id1 = 51919;
id2 = 83340;
targetIDS = [id1, id2];                %Target Stars: single number or vector

% Far Cluster = 50162
% Center Cluster  = 54172

% FASTSHIP FUNCTIONS INPUTS
tm = 1;                    %Lambert input: -1 for long way, +1 for short
t01 = 4.427930866301;
t02 = 8.507563599537;
t0 = [t01, t02];%0:.25:10;             %Launch Time: single number or vector
dmyr1 = 85.572069133699;
dmyr2 = 81.492436400462;
dmyr = [dmyr1,dmyr2]; %[1:.5:80];        %Travel Time: Single number or vector

% RUN MANEUVERS
%[results_table] = fastship_maneuver(starIds, t0Vec, dmyrVec, tm, starData)
[full_results, submission_results] = fastship_maneuver(targetIDS, t0, dmyr, tm, starData);

% MAKE SUBMISSIONS
[qs1, qs2] = fastship_submission_format(submission_results, id1, id2, t01,t02, t01+dmyr1, t02+dmyr2);
submission_output = sprintf('%2.0f,%.0f,%.8f,%.8f,%.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f\n',qs1, qs2)
% fid = fopen('fastship_output.txt','wt');
% fprintf(fid, submission_output);
% fclose(fid);

% MAKE RESULTS TABLE
% [full_results, submission_results] = fastship_maneuver([54172], 0, [8], tm, starData);
%T = array2table(full_results, 'VariableNames',{'LaunchTime','TravelTime','StarID','DV1', 'DV2', 'deltaV'}); %Add r and theta to the table

% MAKE PLOTS
% launchTime = T{:,1};
% %submission_results(submission_results(:,1) == ID1,:);
% travelTime = T{:,2};
% starID = T{:,3};
% totalDV = T{:,6};
% 
% plot3(launchTime,travelTime,totalDV, 'b.')
% xlabel('Launch Time (myr)')
% ylabel('Travel Time (myr)')
% zlabel('Total DV')
% 
% figure
% plot(launchTime, totalDV, 'b.')
% xlabel('Launch Time (myr)')
% ylabel('Total DV')

%function [qs1, qs2] = fastship_submission_format(submission_results, ID1, ID2, t01, t02, dt1, dt2)
[qs1, qs2] = fastship_submission_format(submission_results, 51919, 96507, 7.228994102426, 0, 82+7.228994102426, 80);        
submission_output = sprintf('%2.0f,%.0f,%.10f,%.10f,%.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f\n',qs1, qs2)
fid = fopen('fastship_output.txt','wt');
fprintf(fid, submission_output);
fclose(fid);

%STOP TIMER
toc 


