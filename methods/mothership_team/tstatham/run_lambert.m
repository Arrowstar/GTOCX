clear all
close all
clc

%Start a timer
tic;
 
starData = readmatrix('data/stars.txt');  %Read in Stardata 
ids = starData(:,1);                      %starids
starId0 = 0;

% FASTSHIP FUNCTIONS INPUTS
tm = 1;     %Lambert input: -1 for long way, +1 for short
t0 = 0;     %Launch Time: single number or vector
dmyr = [8, 15]; %Travel Time: Single number or vector
targetIDS = [54172, 50162]; %Target Stars: single number or vector

%[results_table] = fastship_maneuver(starIds, t0Vec, dmyrVec, tm, starData)
[full_results, submission_results] = fastship_maneuver(targetIDS, t0, dmyr, tm, starData);

%[full_results, submission_results] = fastship_maneuver([54172], 0, [8], tm, starData);
T = array2table(full_results, 'VariableNames',{'LaunchTime','TravelTime','StarID','DV1', 'DV2', 'deltaV'}); %Add r and theta to the table

function [qs1, qs2] = fastship_submission_format(submission_results, ID1, ID2, t01, t02, dt1, dt2)

[qs1, qs2] = fastship_submission_format(submission_results, 54172, 50162, 0, 0, 8, 15);
%submission_output = sprintf('%2.0f,%.0f,%.1f,%.1f,%.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f\n',qs1, qs2)

%csvwrite('fastship_submission.csv', qs)
% for i  = 1:1:length(results(:,1))
%     starId = position_table(i,3);
%     rkpc(i) = starData(starId,2);
%     theta(i) = starData(starId,2);
% end


toc %timer


