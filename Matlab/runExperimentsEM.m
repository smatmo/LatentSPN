%
% runExperimentsEM.m
%
% Runs the EM algorithm on 103 datasets and pre-trained SPNs by 
% Poon & Domingos.
%
% We run EM for 
% * any combination of parameters to be trained (weights/means/sigmas)
% * no missing data/33% missing data/66% missing data
% * with random parameter intialization (3 random restart)
% * with original parameters obtained by Poon & Domingos
%
% This code can be executed by several workers, if they operate on the 
% same file system.
%
% Robert Peharz, October 2016
%


clear all
setPaths

binPath = '../CPP/bin/';
dataPath = '../Data/';
modelPath = '../Models/';
resultPath = '../Results/EM/';

caltechModelPath = [PoonDomingosRelPath, '/results/caltech/models/'];
olivettiModelPath = [PoonDomingosRelPath, '/results/olivetti/models/'];

ls = dir(dataPath);
if length(ls) < 3
    error('no datasets found. Did you run convertData.m?')
end

numIter = 30;
INITRAND_range    = [0, 1];
MAKEMISSING_range = [0.0, 0.333, 0.666];
UPDATE_range = [0,0,1;...
    0,1,0;...
    0,1,1;...
    1,0,0;...
    1,0,1;...
    1,1,0;...
    1,1,1];

for irc = 1:length(INITRAND_range)
    INITRAND = INITRAND_range(irc);
    
    for mmc = 1:length(MAKEMISSING_range)
        MAKEMISSING = MAKEMISSING_range(mmc);
        
        for udc = 1:size(UPDATE_range,1)
            UPDATEWEIGHTS = UPDATE_range(udc,3);
            UPDATEMEANS   = UPDATE_range(udc,2);
            UPDATESIGMAS  = UPDATE_range(udc,1);
            
            for k = 3:length(ls)
                dataset = ls(k).name(1:end-4);
                
                fprintf(dataset);
                fprintf('\n')
                clear TestData TestMU TestSIGMA TrainData TrainMU TrainSIGMA
                load([dataPath, dataset, '.mat'])
                
                if strcmp(dataset, 'ORL')
                    inModelFile = [olivettiModelPath, 'olive.mdl'];
                else
                    inModelFile = [caltechModelPath, dataset, '.mdl'];
                end
                
                if MAKEMISSING > 0
                    numD = size(TrainData,2);
                    numMissing = round(numD * MAKEMISSING);
                    
                    if MAKEMISSING > 0.5
                        cf = 10000;
                    else
                        cf = 0;
                    end
                    rng(udc * length(ls) + k + cf);
                    for s = 1:size(TrainData,1)
                        rp = randperm(numD);
                        TrainData(s, rp(1:numMissing)) = NaN;
                    end
                end
                
                if INITRAND
                    for r = 1:3
                        fprintf('restart %d\n',r)
                        
                        if MAKEMISSING == 0
                            resultFile = sprintf('%s/%s_%d_%d_%d_restart%d.mat', resultPath, dataset, UPDATEWEIGHTS, UPDATEMEANS, UPDATESIGMAS, r);
                            outModelFile = sprintf('%s/%s_EM_%d_%d_%d_restart%d.mod',...
                                modelPath, dataset, UPDATEWEIGHTS, UPDATEMEANS, UPDATESIGMAS, r);
                        else
                            resultFile = sprintf('%s/%s_%d_%d_%d_restart%d_missing%f.mat', resultPath, dataset, UPDATEWEIGHTS, UPDATEMEANS, UPDATESIGMAS, r, MAKEMISSING);
                            outModelFile = sprintf('%s/%s_EM_%d_%d_%d_restart%d_missing%f.mod',...
                                modelPath, dataset, UPDATEWEIGHTS, UPDATEMEANS, UPDATESIGMAS, r, MAKEMISSING);
                        end
                        
                        %%% rand seed
                        seed = r * 1000 + udc * 10;
                        
                        if exist(resultFile, 'file')
                            continue;
                        end
                        tmp = 'tmp';
                        save(resultFile, 'tmp', '-v7')
                                                
                        [success, history] = trainSPN_EM(inModelFile, TrainData, TestData, outModelFile, ...
                            'numIter', numIter, ...
                            'updateWeights', UPDATEWEIGHTS, ...
                            'updateMeans', UPDATEMEANS, ...
                            'updateSigmas', UPDATESIGMAS, ...
                            'minSigma', 0.1, ...
                            'stop_relLikelihoodChange', -1e6,...
                            'earlyStoppingK', numIter + 1, ...
                            'PDformat', 1,...
                            'width', 64,...
                            'height', 64,...
                            'randSeed', seed,...
                            'binPath', binPath);                        
                        if success
                            save(resultFile, 'history', '-v7');
                        else
                            delete(resultFile)
                            error('EM training failed')
                        end
                    end
                else
                    if MAKEMISSING == 0
                        resultFile = sprintf('%s/%s_%d_%d_%d.mat', resultPath, dataset, UPDATEWEIGHTS, UPDATEMEANS, UPDATESIGMAS);
                        outModelFile = sprintf('%s/%s_EM_%d_%d_%d.mod',...
                            modelPath, dataset, UPDATEWEIGHTS, UPDATEMEANS, UPDATESIGMAS);
                    else
                        resultFile = sprintf('%s/%s_%d_%d_%d_missing%f.mat', resultPath, dataset, UPDATEWEIGHTS, UPDATEMEANS, UPDATESIGMAS, MAKEMISSING);
                        outModelFile = sprintf('%s/%s_EM_%d_%d_%d_missing%f.mod',...
                            modelPath, dataset, UPDATEWEIGHTS, UPDATEMEANS, UPDATESIGMAS, MAKEMISSING);
                    end
                    
                    if exist(resultFile, 'file')
                        continue;
                    end
                    tmp = 'tmp';
                    save(resultFile, 'tmp', '-v7')
                                        
                    [success, history] = trainSPN_EM(inModelFile, TrainData, TestData, outModelFile, ...
                        'numIter', numIter, ...
                        'updateWeights', UPDATEWEIGHTS, ...
                        'updateMeans', UPDATEMEANS, ...
                        'updateSigmas', UPDATESIGMAS, ...
                        'minSigma', 0.1, ...
                        'stop_relLikelihoodChange', -1e6,...
                        'earlyStoppingK', numIter + 1, ...
                        'PDformat', 1,...
                        'width', 64,...
                        'height', 64,...
                        'binPath', binPath);
                    if success
                        save(resultFile, 'history', '-v7');
                    else                        
                        delete(resultFile)
                        error('EM training failed')
                    end
                end
            end
        end
    end
end
