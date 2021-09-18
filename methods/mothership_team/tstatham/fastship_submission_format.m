function [qs1, qs2] = fastship_submission_format(submission_results, ID1, ID2, t01, t02, dt1, dt2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
% submission_results = target stars for the lambert solver. Should match stars.txt file
% ID1/ID2 = starID of the two fastship targets
% t01/t02 = launch time of the two fastships
% dt1/dt2 = travel time of the two fastships
%
%OUTPUTS:
% Creates CSV file for submission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%format bank
format long g

%First Fastship result
submission1 = submission_results(submission_results(:,1) == ID1,:);
submission1 = submission1(submission1(:,2) == t01,:);
submission1 = submission1(submission1(:,3) == dt1,:);


%Second Fastship Result
submission2 = submission_results(submission_results(:,1) == ID2,:);
submission2 = submission2(submission2(:,2) == t02,:);
submission2 = submission2(submission2(:,3) == dt2,:);

qs1 = [-11, submission1];
qs2 = [-12, submission2];

submission_output = sprintf('%2.0f,%.0f,%.1f,%.1f,%.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f\n',qs1, qs2);
