function [success, history] = trainSPN_EM(modelFile, Data, valData, outputFile, varargin)
%
% This is a Matlab wrapper for calling the binary EM training tool.
%

success = false;

%%%%%%%%%%%%%%
%%% params %%%
%%%%%%%%%%%%%%

numIter           = [];
stop_relLikelihoodChange = [];
earlyStoppingK    = [];
updateWeights     = [];
updateMeans       = [];
updateSigmas      = [];
minSigma          = [];
randSeed          = [];
trainDataFile     = [];
valDataFile       = [];
outputLabels      = [];
binPath           = '';
PDformat          = [];
width             = [];
height            = [];

k = 1;
while k < length(varargin)
    switch varargin{k}
        case 'numIter'
            numIter = varargin{k+1};
            k = k+1;
        case 'stop_relLikelihoodChange'
            stop_relLikelihoodChange = varargin{k+1};
            k = k+1;
         case 'earlyStoppingK'
            earlyStoppingK = varargin{k+1};
            k = k+1;
        case 'updateWeights'
            updateWeights = varargin{k+1};
            k = k+1;
        case 'updateMeans'
            updateMeans = varargin{k+1};
            k = k+1;
        case 'updateSigmas'
            updateSigmas = varargin{k+1};
            k = k+1;
        case 'minSigma'
            minSigma = varargin{k+1};
            k = k+1;
        case 'randSeed'
            randSeed = varargin{k+1};
            k = k+1;
        case 'trainDataFile'
            trainDataFile = varargin{k+1};
            k = k+1;
        case 'valDataFile'
            valDataFile = varargin{k+1};
            k = k+1;
        case 'binPath'
            binPath = varargin{k+1};
            k = k+1;
        case 'PDformat'
            PDformat = varargin{k+1};
            k = k+1;
        case 'width'
            width = varargin{k+1};
            k = k+1;
        case 'height'
            height = varargin{k+1};
            k = k+1;
        otherwise
            disp(varargin{k})
            error('unknown parameter')
    end
    k = k + 1;
end

if ~exist(sprintf('%strainEM', binPath), 'file')
    error('binary trainEM not found. Did you compile the C++ source?')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(valData) && isempty(valDataFile)
    if size(Data,2) ~= size(valData,2)
        error('size(Data,2) ~= size(valData,2)')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpFile = tempname;

if isempty(trainDataFile)
    tmpTrainFile = [tmpFile, '_trainSPN_EM'];
    writeSPNVecList(tmpTrainFile, Data);
else
    tmpTrainFile = trainDataFile;
end

cmd = sprintf('"%strainEM" %s %s %s ',...
    binPath, modelFile, tmpTrainFile, outputFile);

if isempty(valDataFile)
    if ~isempty(valData)
        tmpValFile = [tmpFile, '_valSPN_PD'];
        writeSPNVecList(tmpValFile, valData);
        cmd = [cmd, sprintf(' valDataFile = %s', tmpValFile)];
    end
else
    tmpValFile = valDataFile;
    cmd = [cmd, sprintf(' valDataFile = %s', tmpValFile)];
end

if ~isempty(numIter)
    cmd = [cmd, sprintf(' numIter = %d', numIter)];
end

if ~isempty(stop_relLikelihoodChange)
    cmd = [cmd, sprintf(' stop_relLikelihoodChange = %d', stop_relLikelihoodChange)];
end

if ~isempty(earlyStoppingK)
    cmd = [cmd, sprintf(' earlyStoppingK = %d', earlyStoppingK)];
end

if ~isempty(updateWeights)
    cmd = [cmd, sprintf(' updateWeights = %d', updateWeights)];
end

if ~isempty(updateMeans)
    cmd = [cmd, sprintf(' updateMeans = %d', updateMeans)];
end

if ~isempty(updateSigmas)
    cmd = [cmd, sprintf(' updateSigmas = %d', updateSigmas)];
end

if ~isempty(minSigma)
    cmd = [cmd, sprintf(' minSigma = %d', minSigma)];
end

if ~isempty(randSeed)
    cmd = [cmd, sprintf(' randSeed = %d', randSeed)];
end

if ~isempty(PDformat)
    cmd = [cmd, sprintf(' PDformat = %d', PDformat)];
end

if ~isempty(width)
    cmd = [cmd, sprintf(' width = %d', width)];
end

if ~isempty(height)
    cmd = [cmd, sprintf(' height = %d', height)];
end

if nargout > 1
    historyFile = [tmpFile, '_hist'];
    cmd = [cmd, sprintf(' historyFile = %s', historyFile)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sysStat = system(cmd);

%%% clean up
if isempty(trainDataFile)
    delete(tmpTrainFile);
end

if isempty(valDataFile) && ~isempty(valData)
    delete(tmpValFile);
end

if ~isempty(outputLabels)
    delete(tmpLabelFile);
end

if (sysStat ~= 0)
    history = [];
    success = false;
    return;
end

if nargout > 1
    history = readSPNVecList(historyFile);
    delete(historyFile);
end

success = true;

