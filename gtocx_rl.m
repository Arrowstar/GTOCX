clc; clear all; format long g; close all;

addpath(genpath('.'));
warning('off');

env = GTOCXEnvironment('stars.txt');
obsInfo = getObservationInfo(env);
actInfo = getActionInfo(env);

validateEnvironment(env);

agent = rlTD3Agent(obsInfo, actInfo);
% agent.AgentOptions.SampleTime = 0.1;
trainOpts = rlTrainingOptions('MaxEpisodes',100000, 'MaxStepsPerEpisode',1000, 'ScoreAveragingWindowLength',10, 'StopTrainingCriteria','AverageReward', ...
                              'StopTrainingValue',1000, 'UseParallel',false, "Plots","none");
trainOpts.Verbose = true;
% trainOpts.ParallelizationOptions.Mode = "async";

train(agent, env, trainOpts);