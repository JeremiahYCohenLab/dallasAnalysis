function createStartValues(filename, paramRanges)

if nargin < 2 
    paramRanges = [0.1 1; 0.1 1; 0.1 1; 1 9; 0 1];
end
    
for i = 1:length(paramRanges)
    temp = linspace(paramRanges(i,1), paramRanges(i,2), 100);
    temp = repmat(temp, [1,5]);
    temp = temp(randperm(length(temp)));
    startVals(:,i) = temp';
end

if strfind(filename, 'opponency')
    savepath = strcat('C:\Users\cooper_PC\Desktop\githubRepositories\cooperAnalysis\matlabCode\operantMatching\learningModels\qLearning\opponency\', filename);
elseif strfind(filename, 'dual')
    savepath = strcat('C:\Users\cooper_PC\Desktop\githubRepositories\cooperAnalysis\matlabCode\operantMatching\learningModels\qLearning\dual\', filename);
elseif strfind(filename, 'spike')
    savepath = strcat('C:\Users\cooper_PC\Desktop\githubRepositories\cooperAnalysis\matlabCode\operantMatching\functions\spikeAnalysis\', filename);
elseif strfind(filename, 'LR')
    savepath = strcat('C:\Users\cooper_PC\Desktop\githubRepositories\cooperAnalysis\matlabCode\operantMatching\learningModels\qLearning\avgLearningRate\', filename);
else
    savepath = strcat('C:\Users\cooper_PC\Desktop\githubRepositories\cooperAnalysis\matlabCode\operantMatching\learningModels\qLearning\', filename);
end

csvwrite([savepath], startVals)