%
% plotLLs.m
%
% Plots the likelihoods of EM 
%
% Robert Peharz, October 2016
%

clear all
close all

resultPath = '../Results/EM/';
dataPath = '../Data/';

% which param combinations shall be plotted?
% [sigmas, means, weights]
updateFlags = ...
    [
    [0,0,1]
    [0,1,0]
    [0,1,1]
    [1,0,0]
    [1,0,1]
    [1,1,0]
    [1,1,1]    
    ];

%%% the following settings reproduce Fig. 9 in the paper
% includeNumRandomInits = 3;
% includeOriginalInits = 0; 
% includeMAKEMISSING = [0.0];

includeNumRandomInits = 3;  % 0-3: number of random restarts to include 
includeOriginalInits = 0;   % 1 to include Poon & Domingos original params
includeMAKEMISSING = [0.0]; % includeMAKEMISSING = [0.0, 0.333, 0.666];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TrainLLs = cell(size(updateFlags,1),1);
TestLLs = cell(size(updateFlags,1),1);


%%% get data set sizes
ls = dir(dataPath);
DATASETS = cell(length(ls)-2, 1);
TRAINDATA_N = zeros(length(ls)-2, 1);
TESTDATA_N = zeros(length(ls)-2, 1);

for k=3:length(ls)
    load([dataPath, '/', ls(k).name]);
    trainN = size(TrainData,1);
    testN = size(TestData,1);
    DATASETS{k-2} = ls(k).name(1:end-4);
    TRAINDATA_N(k-2) = trainN;
    TESTDATA_N(k-2) = testN;
end


%%% load results
ls = dir(resultPath);
numResUsed = 0;
for k=3:length(ls)
    
    fileName = ls(k).name;
    
    %%% filter MISSING
    idx = strfind(fileName, 'missing');
    if idx
        MISSING = str2double(fileName((idx(end)+7):end-4));
        fileName = fileName(1:idx(end)-2);
    else
        MISSING = 0.0;
        fileName = fileName(1:end-4);
    end
    
    if ~any(MISSING == includeMAKEMISSING)
        continue
    end
    
    
    %%% filter inits
    idx = strfind(fileName, 'restart');
    if idx
        restartNum = str2double(fileName(idx(end)+7));
        fileName = fileName(1:idx(end)-2);
        if restartNum > includeNumRandomInits
            continue
        end
    else
        restartNum = -1;
        if includeOriginalInits == 0
            continue
        end
    end
    
   
    %%% filter param switch
    updateSigmas = str2double(fileName(end));
    fileName = fileName(1:end-2);
    updateMeans = str2double(fileName(end));
    fileName = fileName(1:end-2);
    updateWeights = str2double(fileName(end));
    fileName = fileName(1:end-2);
    
    curPS = [updateSigmas,updateMeans,updateWeights];
    
    uf = all(updateFlags == repmat(curPS, size(updateFlags,1),1),2);
    if ~any(uf)
        continue
    end
    ufIdx = find(uf,1);    
    
    
    %%% dataset name
    dataset = fileName;
       
    clear tmp history
    load([resultPath, ls(k).name])
    if exist('tmp', 'var')
        warning('tmp file found')
        continue
    end
    
    trainN = -1;
    testN = -1;
    for l=1:length(DATASETS)
        if strcmp(DATASETS(l), dataset)
            trainN = TRAINDATA_N(l);
            testN = TESTDATA_N(l);
            break
        end
    end
    if (trainN < 0) || (testN < 0)
        error('no dataset size found');
    end
    
    numResUsed = numResUsed + 1;
    
    TrainLLs{ufIdx} = [TrainLLs{ufIdx}; history{2} / trainN];
    TestLLs{ufIdx} = [TestLLs{ufIdx}; history{3} / testN];
end

numResExp = size(updateFlags,1) * ...
    (includeNumRandomInits + (includeOriginalInits~=0)) * ...
    length(includeMAKEMISSING) * 103;

if numResExp ~= numResUsed
    warning('expected number of results does not match with used number of results')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% plotting styles
cols = {'b','r','g',[0.9,0.45,0],'k','m','c'};
styls = {'*','s','p','v','+','<','>'};

%%% make texts for the legen
legendTexts = cell(size(updateFlags,1),1);
for k=1:size(updateFlags,1)
    t = '';
    if updateFlags(k,3)
        t='W';
    end
    if updateFlags(k,2)
        t=[t,'M'];
    end
    if updateFlags(k,1)
        t=[t,'V'];
    end
    
    legendText{k} = t;    
end

%%%%%%%%%%%%%%%%%%%%
%%% training set %%%
%%%%%%%%%%%%%%%%%%%%

figure(1);
clf
hold on

%%% plot markers
for f = 1:length(TrainLLs)   
    plot(0:4:30,mean(TrainLLs{f}(:,1:4:end)), styls{f}, 'linewidth', 1.5, 'Color', cols{f})
end

legend(legendText);

%%% plot lines
for f = 1:length(TrainLLs)   
    plot(0:30,mean(TrainLLs{f}), 'linewidth', 1.5, 'Color', cols{f})
end

xlabel('Iteration')
ylabel('Train log-likelihood')
grid on
box on


%%%%%%%%%%%%%%%%
%%% test set %%%
%%%%%%%%%%%%%%%%

figure(2);
clf
hold on

%%% plot markers
for f = 1:length(TrainLLs)   
    plot(0:4:30,mean(TestLLs{f}(:,1:4:end)), styls{f}, 'linewidth', 1.5, 'Color', cols{f})
end

legend(legendText);

%%% plot lines
for f = 1:length(TestLLs)   
    plot(0:30,mean(TestLLs{f}), 'linewidth', 1.5, 'Color', cols{f})
end

xlabel('Iteration')
ylabel('Test log-likelihood')
grid on
box on

