clc; clear all; format long g; close all;

addpath(genpath('.'));

env = GTOCXEnvironment('stars.txt');
obsInfo = getObservationInfo(env);
actInfo = getActionInfo(env);

validateEnvironment(env);