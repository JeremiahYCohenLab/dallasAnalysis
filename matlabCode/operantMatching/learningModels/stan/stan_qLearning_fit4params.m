function [fit, samples] = stan_qLearning_fit4params(xlFile, animal, category, revForFlag)

if nargin < 4
    revForFlag = 0;
end

[root, sep] = currComputer();

[~, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

session = zeros(630, length(dayList));
trial = zeros(630, length(dayList));
choice = zeros(630, length(dayList));
outcome = zeros(630, length(dayList));

for i = 1:length(dayList)
    sessionName = dayList{i};
    filename = [sessionName '.asc'];
    if revForFlag == 1
        [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData_revFor(filename);
        behavStruct = parseBehavioralData_revFor(behSessionData, unCorrectedBlockSwitch);
    else
        [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename);
        behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);
    end

    choiceTmp = behavStruct.allChoices;
    choiceTmp(choiceTmp == -1) = 2;
    choiceTmp = choiceTmp';
    outcomeTmp = abs(behavStruct.allRewards); 
    outcomeTmp = outcomeTmp';
    sessionTmp = repmat(i, length(outcomeTmp), 1);
    trialTmp = [1:length(outcomeTmp)]';
    
    session(1:trialTmp(end),i) = sessionTmp;
    trial(1:trialTmp(end),i) = trialTmp;
    choice(1:trialTmp(end),i) = choiceTmp;
    outcome(1:trialTmp(end),i) = outcomeTmp;
    Tsesh(i,1) = trialTmp(end);
    
end

N = length(dayList);
T = 630;


session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'choice', choice', 'outcome', outcome');

fit = stan('file','stan_qLearning_4params.stan','data',session_dat,'verbose',true);

print(fit);

samples = fit.extract('permuted',true);
fprintf('\n\tmean    sd\n')
fprintf('aN\t%1.3f\t%1.3f\n',mean(samples.aN),std(samples.aN));
fprintf('aP\t%1.3f\t%1.3f\n',mean(samples.aP),std(samples.aP));
fprintf('beta\t%1.3f\t%1.3f\n',mean(samples.beta),std(samples.beta));
fprintf('mu_aN\t%1.3f\t%1.3f\n',mean(samples.mu_aN),std(samples.mu_aN));
fprintf('mu_aP\t%1.3f\t%1.3f\n',mean(samples.mu_aP),std(samples.mu_aP));
fprintf('mu_beta\t%1.3f\t%1.3f\n',mean(samples.mu_beta),std(samples.mu_beta));

