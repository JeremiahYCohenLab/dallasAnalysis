% The analyses below reproduce Fig. 2b in:
%   Cohen JY, Haesler S, Vong L, Lowell BB, Uchida N. Neuron-type-specific
%     signals for reward and punishment in the ventral tegmental area.
%     Nature 482: 85-88, 2012.
% The following R code saves the file brAll.rda as a Matlab data file,
%   using the library (R.matlab); the advantage is that the R data file
%   is only about 3 Mb:
%   writeMat('brAll.mat', brAll)
% Written by Jeremiah Y. Cohen


% Load data
% This is a 95 x 234 x 559 array, corresponding to 95 neurons, a maximum of
%   234 trials, and a maximum of 559 spike times per trial
% Spike times are relative to CS onset, within a window of 1 s before CS
%   onset to 4 s after, expressed in ms
% The first 49 neurons are of type I, the next 29 are of type II, the last
%   17 are of type III
load('brAll.mat');

% Generate spike density functions for each neuron
mySDFs = [];
h = waitbar(0,'Please wait...');
k = 0;
for i = 1 : 4: 95
    for j = 1 : 234
        temp = brAll(i, j, :);
        if sum(~isnan(temp)) > 0
            mySDFs = cat(3, mySDFs, sdf(timesToCounts(temp)));
        end
        k = k+1;
        waitbar(k/(24*234))
    end
end

% Calculate auROC for each neuron for big-reward trials
auROCdata = [];
for i = 1 : 95
    auROCdata = cat(2, auROCdata, auROC(mySDFs(i)));
end

% Calculate PCA on auROC
myPCA = princomp(auROCdata);

% Cluster neurons based on PCA
myCluster = linkage(myPCA, 'complete', 'euclidean');
