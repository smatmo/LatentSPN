%
% deleteDummyFiles.m
%
% Delete placeholder files used for parallel EM training.
%
% If workers crash during EM training, they might leave "dummy" results
% files in the results folder. This script deletes these dummy files.
%
% Robert Peharz, October 2016
%

clear all
close all

resultPath = '../Results/EM/';

ls = dir(resultPath);
numResUsed = 0;
for k=3:length(ls)
    clear tmp history
    load([resultPath, ls(k).name])
    if exist('tmp', 'var')
        fprintf('delete %s\n', ls(k).name);
        delete([resultPath, ls(k).name]);
    end    
end
