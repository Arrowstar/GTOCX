clc; clear all; format long g; close all;

addpath(genpath('.'));

useParallel = false;

if(useParallel && isempty(gcp('nocreate')))
    parpool(feature('numcores'));
    parfevalOnAll(gcp(), @warning, 0, 'off');
end

env = GTOCXEnvironment('stars.txt');
obsInfo = getObservationInfo(env);
actInfo = getActionInfo(env);

validateEnvironment(env);

agent = rlTD3Agent(obsInfo, actInfo);
agent.AgentOptions.ExplorationModel.LowerLimit = -1;
agent.AgentOptions.ExplorationModel.UpperLimit = 1;

trainOpts = rlTrainingOptions('MaxEpisodes',100000, 'MaxStepsPerEpisode',1000, 'ScoreAveragingWindowLength',10, 'StopTrainingCriteria','GlobalStepCount', ...
                              'StopTrainingValue',1E99, 'UseParallel',useParallel, "Plots","none");
trainOpts.Verbose = true;
% trainOpts.ParallelizationOptions.Mode = "async";

train(agent, env, trainOpts);