%
% convertData.m
%
% This file loads the data and processes it the same way as by 
% Poon & Domingos, 2011:
% - take the last third (but at most 50) samples as test data
% - normalize zero mean/unit variance per sample
% - take the central 64x64 patch of the Caltech images
% Then the data is stored in Matlab format.
%
% Robert Peharz, October 2016
%

clear all
setPaths


%%%
%%% Convert Caltech
%%%
CaltechBasePath = [PoonDomingosRelPath, 'data/caltech/'];
if ~exist(CaltechBasePath, 'dir')
    error('Caltech path not found. Have you downloaded and added the Poon & Domingos package?')
end
ls = dir(CaltechBasePath);

for k = 3:length(ls)
    if ~ls(k).isdir
        continue
    end
    
    lss = dir([CaltechBasePath, ls(k).name]);
    
    Data = [];
    MU = [];
    SIGMA = [];
    for l=3:length(lss)
        im = dlmread([CaltechBasePath, ls(k).name,'/',lss(l).name]);
        im = im(19:82,19:82);
                
        mu = mean(im(:));
        sigma = std(im(:));
        im = (im - mu) / sigma;
        
        Data = [Data; im(:)'];
        MU = [MU; mu];
        SIGMA = [SIGMA; sigma];
    end
   
    N = size(Data,1);
    testN = min(50, floor(N / 3));
    trainN = N - testN;
    
    TrainData = Data(1:trainN,:);
    TrainMU = MU(1:trainN);
    TrainSIGMA = SIGMA(1:trainN);
    TestData = Data(trainN+1:end,:);
    TestMU = MU(trainN+1:end);
    TestSIGMA = SIGMA(trainN+1:end);
    
    save(['../Data/', ls(k).name, '.mat'], 'TrainData', 'TestData', 'TrainMU', 'TrainSIGMA', 'TestMU', 'TestSIGMA');
end



%%%
%%% Convert Olivetti
%%%
ORLPath = [PoonDomingosRelPath, 'data/olivetti/olivetti.raw'];
if ~exist(ORLPath, 'file')
    error('ORL file not found. Have you downloaded and added the Poon & Domingos package?')
end
DataRaw = dlmread(ORLPath);

Data = [];
MU = [];
SIGMA = [];
for k=1:size(DataRaw,2)
    im = reshape(DataRaw(:,k), 64, 64);
    
    mu = mean(im(:));
    sigma = std(im(:));
    im = (im - mu) / sigma;
    
    Data = [Data; im(:)'];
    MU = [MU; mu];
    SIGMA = [SIGMA; sigma];    
end

N = size(Data,1);
testN = min(50, floor(N / 3));
trainN = N - testN;

TrainData = Data(1:trainN,:);
TrainMU = MU(1:trainN);
TrainSIGMA = SIGMA(1:trainN);
TestData = Data(trainN+1:end,:);
TestMU = MU(trainN+1:end);
TestSIGMA = SIGMA(trainN+1:end);

save('../Data/ORL.mat', 'TrainData', 'TestData', 'TrainMU', 'TrainSIGMA', 'TestMU', 'TestSIGMA');

