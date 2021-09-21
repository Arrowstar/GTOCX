clc; clear all; format long g; close all;

addpath(genpath('.'));

env = GTOCXEnvironment('stars.txt');
obsInfo = getObservationInfo(env);
actInfo = getActionInfo(env);

validateEnvironment(env);

agent = rlPPOAgent(obsInfo, actInfo);
% agent.AgentOptions.SampleTime = 0.1;
trainOpts = rlTrainingOptions('MaxEpisodes',100000, 'MaxStepsPerEpisode',100, 'ScoreAveragingWindowLength',10, 'StopTrainingCriteria','GlobalStepCount', ...
                              'StopTrainingValue',1E99, 'UseParallel',false, "Plots","none");
trainOpts.Verbose = true;
% trainOpts.ParallelizationOptions.Mode = "async";

train(agent, env, trainOpts);