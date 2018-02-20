%
% runExperimentsMPE.m
%
% Runs MPE inference and displays results.
% The actual experiments are done in the C++ binary MPEsynth.
%
% Robert Peharz, October 2016
%

clear all

DRange = [2,3,4];
alphaRange = [0.5,1,2.0];

binPath = '../CPP/bin/';
resultPath = '../Results/MPE/';

%%% comment out the following lines, 
%%% if you don't want to re-run the MPE experiment
if ~exist(sprintf('%sMPEsynth', binPath), 'file')
    error('binary MPEsynth not found. Did you compile the C++ source?')
end
cmd = sprintf('%s/MPEsynth %s', binPath, resultPath);
sysStat = system(cmd);
%%%
%%%

for DRangeCount = 1:length(DRange)
    D = DRange(DRangeCount);
    for aRangeCount = 1:length(alphaRange)
        alpha = alphaRange(aRangeCount);
        
        filename = sprintf('%s/MPEresult_D%d_alpha%d.txt', resultPath, DRangeCount-1, aRangeCount-1);
        R = dlmread(filename);
        
        LLorig{DRangeCount,aRangeCount} = R(:,1:3);
        LLaug_uni{DRangeCount,aRangeCount} = R(:,4:6);
        LLaug_det{DRangeCount,aRangeCount} = R(:,7:9);
        
        
        %%% sanity checks
        if any(LLorig{DRangeCount,aRangeCount}(:,1) < LLorig{DRangeCount,aRangeCount}(:,2) - 1e-12)
            error('exhaustive solution sub-optimal')
        end
        if any(LLorig{DRangeCount,aRangeCount}(:,1) < LLorig{DRangeCount,aRangeCount}(:,3) - 1e-12)
            error('exhaustive solution sub-optimal')
        end
        
        if any(LLaug_uni{DRangeCount,aRangeCount}(:,1) < LLaug_uni{DRangeCount,aRangeCount}(:,2) - 1e-12)
            error('exhaustive solution sub-optimal')
        end
        if any(LLaug_uni{DRangeCount,aRangeCount}(:,1) < LLaug_uni{DRangeCount,aRangeCount}(:,3) - 1e-12)
            error('exhaustive solution sub-optimal')
        end
        
        if any(LLaug_det{DRangeCount,aRangeCount}(:,1) < LLaug_det{DRangeCount,aRangeCount}(:,2) - 1e-12)
            error('exhaustive solution sub-optimal')
        end
        if any(LLaug_det{DRangeCount,aRangeCount}(:,1) < LLaug_det{DRangeCount,aRangeCount}(:,3) - 1e-12)
            error('exhaustive solution sub-optimal')
        end
        %%%
        
    end
end


fprintf('\n');
fprintf('In augmented SPN, uniform weights\n');
fprintf('--------------------------------------------------------\n');
fprintf('                       UNI          DET\n');
fprintf('--------------------------------------------------------\n');
for DRangeCount = 1:length(DRange)
    D = DRange(DRangeCount);
    for aRangeCount = 1:length(alphaRange)
        alpha = alphaRange(aRangeCount);
        
        LLtmp = LLaug_uni{DRangeCount, aRangeCount};
        numCorrectUni = sprintf('%d', sum(LLtmp(:,2) > LLtmp(:,1) - 1e-12));
        numCorrectUni = [repmat(' ',1,3-length(numCorrectUni)), numCorrectUni];
        numCorrectDet = sprintf('%d', sum(LLtmp(:,3) > LLtmp(:,1) - 1e-12));
        numCorrectDet = [repmat(' ',1,3-length(numCorrectDet)), numCorrectDet];
        meanUni = mean(LLtmp(:,2) - LLtmp(:,1));
        meanDet = mean(LLtmp(:,3) - LLtmp(:,1));
        
        fprintf('D = %d   alpha = %0.2f   %0.2f (%s)   %0.2f (%s)\n', D, alpha, meanUni, numCorrectUni, meanDet, numCorrectDet);
    end
    fprintf('--------------------------------------------------------\n');
end



fprintf('\n');
fprintf('In augmented SPN, deterministic weights\n');
fprintf('--------------------------------------------------------\n');
fprintf('                       UNI          DET\n');
fprintf('--------------------------------------------------------\n');
for DRangeCount = 1:length(DRange)
    D = DRange(DRangeCount);
    for aRangeCount = 1:length(alphaRange)
        alpha = alphaRange(aRangeCount);
        
        LLtmp = LLaug_det{DRangeCount, aRangeCount};
        numCorrectUni = sprintf('%d', sum(LLtmp(:,2) > LLtmp(:,1) - 1e-12));
        numCorrectUni = [repmat(' ',1,3-length(numCorrectUni)), numCorrectUni];
        numCorrectDet = sprintf('%d', sum(LLtmp(:,3) > LLtmp(:,1) - 1e-12));
        numCorrectDet = [repmat(' ',1,3-length(numCorrectDet)), numCorrectDet];
        meanUni = mean(LLtmp(:,2) - LLtmp(:,1));
        meanDet = mean(LLtmp(:,3) - LLtmp(:,1));
        
        fprintf('D = %d   alpha = %0.2f   %0.2f (%s)   %0.2f (%s)\n', D, alpha, meanUni, numCorrectUni, meanDet, numCorrectDet);
    end
    fprintf('--------------------------------------------------------\n');
end



fprintf('\n');
fprintf('In original SPN\n');
fprintf('--------------------------------------------------------\n');
fprintf('                       UNI          DET\n');
fprintf('--------------------------------------------------------\n');
for DRangeCount = 1:length(DRange)
    D = DRange(DRangeCount);
    for aRangeCount = 1:length(alphaRange)
        alpha = alphaRange(aRangeCount);
        
        LLtmp = LLorig{DRangeCount, aRangeCount};
        numCorrectUni = sprintf('%d', sum(LLtmp(:,2) > LLtmp(:,1) - 1e-12));
        numCorrectUni = [repmat(' ',1,3-length(numCorrectUni)), numCorrectUni];
        numCorrectDet = sprintf('%d', sum(LLtmp(:,3) > LLtmp(:,1) - 1e-12));
        numCorrectDet = [repmat(' ',1,3-length(numCorrectDet)), numCorrectDet];
        meanUni = mean(LLtmp(:,2) - LLtmp(:,1));
        meanDet = mean(LLtmp(:,3) - LLtmp(:,1));
        
        fprintf('D = %d   alpha = %0.2f   %0.2f (%s)   %0.2f (%s)\n', D, alpha, meanUni, numCorrectUni, meanDet, numCorrectDet);
    end
    fprintf('--------------------------------------------------------\n');
end


