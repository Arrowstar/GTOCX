clc; clear all; format long g; close all;

addpath(genpath('.'));
warning('off');

env = GTOCXEnvironment('stars.txt');
obsInfo = getObservationInfo(env);
actInfo = getActionInfo(env);

validateEnvironment(env);

agent = rlPPOAgent(obsInfo, actInfo);
trainOpts = rlTrainingOptions('MaxEpisodes',10000, 'MaxStepsPerEpisode',100, 'ScoreAveragingWindowLength',10, 'StopTrainingCriteria','AverageReward', ...
                              'StopTrainingValue',1000, 'UseParallel',false, "Plots","none");
trainOpts.Verbose = true;
% trainOpts.ParallelizationOptions.Mode = "async";

train(agent, env, trainOpts);