%
% checkMonotonicity.m
%
% Check if EM shows monotonicity, for all EM runs.
% When EM has converged and keeps the likelihood actually unchanged,
% we sometimes observe tiny decreases due to numical artifacts.
% If such a case occurs, we plot the training curve for 1 second,
% indicating decreases as red 'x'.
%
% Robert Peharz, October 2016
%

path = '../Results/EM/';
ls = dir(path);

minD = inf;

for k=3:length(ls)
    
    disp(ls(k).name)
    
    clear history
    load([path, ls(k).name]);
    
    trainLL = history{2};
    D = trainLL(2:end) - trainLL(1:end-1);
    minD = min(minD, min(D));
    
    idx = find(D < 0);
    if ~isempty(idx)
        figure(1)
        clf
        hold on
        plot(trainLL)
        plot(idx+1, trainLL(idx+1), 'rx');
        drawnow
        tic
        while toc < 1
        end
    end
end

fprintf('\nMinimal increase: %d\n\n', minD);

