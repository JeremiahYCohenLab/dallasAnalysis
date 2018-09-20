function models = modelFits_operantMatching(filename, testFlag, saveFigFlag)
if nargin < 2
    testFlag = 0;
    saveFigFlag = 1;
end
if nargin < 3
    saveFigFlag = 1;
end

[sessionData, unCorrectedBlockSwitch, ~, out] = loadBehavioralData(filename);
behavStruct = parseBehavioralData(sessionData, unCorrectedBlockSwitch);
[root, sep] = currComputer_operantMatching();

outcome = abs([behavStruct.allReward_R; behavStruct.allReward_L])';
choice = abs([behavStruct.allChoice_R; behavStruct.allChoice_L])';

% Fit all data
[qFit] = qLearning_fitAll(choice, outcome, testFlag); sprintf('Q-learning finished')
[daFit]= directActor_fitAll(choice, outcome, testFlag); sprintf('Direct actor finished')
[davFit] = directActorWithValue_fitAll(choice, outcome, testFlag); sprintf('Direct actor with value finished')
[acfFit] = actorCritic_fukai_fitAll(choice, outcome, testFlag); sprintf('Actor-critic Fukai finished')
[dtsFit] = dynamicThompsonSampling_fitAll(choice, outcome, testFlag, 1); sprintf('Dynamic Thompson sampling finished')
[spFit] = stackedProbability_fitAll(choice, outcome, testFlag); sprintf('Stacked probability finished')
[phFit] = pearceHall_fitAll(choice, outcome, testFlag); sprintf('Pearce-Hall finished')

models.qLearning = qFit;
models.directActor = daFit;
models.directActorWithValue = davFit;
models.actorCriticFukai = acfFit;
models.dynamicThompson = dtsFit;
models.stackedProb = spFit;
models.pearceHall = phFit;

%% Determine best model within each model class
% Q learning
qLearningModels = fields(qFit);
qBIC = NaN(1,length(qLearningModels));
qLH = NaN(1, length(qLearningModels));
for i = 1:length(qLearningModels)
    qBIC(i) = qFit.(qLearningModels{i}).BIC;
    qLH(i) = qFit.(qLearningModels{i}).LH;
end
[~, best_qLearningModel] = min(qBIC);
probChoice_qLearning = qFit.(qLearningModels{best_qLearningModel}).probChoice;

% Direct actor
directActorModels = fields(daFit);
daBIC = NaN(1,length(directActorModels));
daLH = NaN(1, length(directActorModels));
for i = 1:length(directActorModels)
    daBIC(i) = daFit.(directActorModels{i}).BIC;
    daLH(i) = daFit.(directActorModels{i}).LH;
end
[~, best_directActorModel] = min(daBIC);
probChoice_directActor = daFit.(directActorModels{best_directActorModel}).probChoice;

% Direct actor with value
directActorWithValueModels = fields(davFit);
davBIC = NaN(1, length(directActorWithValueModels));
davLH = NaN(1, length(directActorWithValueModels));
for i = 1:length(directActorWithValueModels)
    davBIC(i) = davFit.(directActorWithValueModels{i}).BIC;
    davLH(i) = davFit.(directActorWithValueModels{i}).LH;
end
[~, best_directActorWithValueModel] = min(davBIC);
probChoice_directActorWithValue = davFit.(directActorWithValueModels{best_directActorWithValueModel}).probChoice;

% Actor critic fukai
acfModels = fields(acfFit);
acfBIC = NaN(1, length(acfModels));
acfLH = NaN(1, length(acfModels));
for i = 1:length(acfModels)
    acfBIC(i) = acfFit.(acfModels{i}).BIC;
    acfLH(i) = acfFit.(acfModels{i}).LH;
end
[~, best_acfModel] = min(acfBIC);
probChoice_acf = acfFit.(acfModels{best_acfModel}).probChoice;

% Dynamic Thompson sampling
dynamicThompsonModels = fields(dtsFit);
dtsBIC = NaN(1, length(dynamicThompsonModels));
dtsLH = NaN(1, length(dynamicThompsonModels));
for i = 1:length(dynamicThompsonModels)
    dtsBIC(i) = dtsFit.(dynamicThompsonModels{i}).BIC;
    dtsLH(i) = dtsFit.(dynamicThompsonModels{i}).LH;
end
[~, best_dynamicThompsonModel] = min(dtsBIC);
probChoice_dynamicThompson = dtsFit.(dynamicThompsonModels{best_dynamicThompsonModel}).probChoice;

% Stacked probability model
stackedProbModels = fields(spFit);
spBIC = NaN(1, length(stackedProbModels));
spLH = NaN(1, length(stackedProbModels));
for i = 1:length(stackedProbModels)
    spBIC(i) = spFit.(stackedProbModels{i}).BIC;
    spLH(i) = spFit.(stackedProbModels{i}).LH;
end
[~, best_stackedProbModel] = min(spBIC);
probChoice_stackedProb = spFit.(stackedProbModels{best_stackedProbModel}).probChoice;

% Pearce-Hall model
phModels = fields(phFit);
phBIC = NaN(1, length(phModels));
phLH = NaN(1, length(phModels));
for i = 1:length(phModels)
    phBIC(i) = phFit.(phModels{i}).BIC;
    phLH(i) = phFit.(phModels{i}).LH;
end
[~, best_phModel] = min(phBIC);
probChoice_ph = phFit.(phModels{best_phModel}).probChoice;

%% Plot

figure; suptitle(filename)
fig_rawData = subplot(5,1,1);
fig_probAction = subplot(5,1,2);
fig_compareModelAndSmoothBehavior = subplot(5, 1, 3);
fig_LH = subplot(5, 1, 4);
fig_BICscores = subplot(5, 1, 5);

xvals = behavStruct.responseInds;

colors = linspecer(6,'qualitative');
qLearning_color = colors(1,:);
directActor_color = colors(2,:);
acf_color = colors(3,:);
dynamicThompson_color = colors(4,:);
stackedProb_color = colors(5,:);
ph_color = colors(6, :);

tMax = length(sessionData);

% Raw data plot
subplot(fig_rawData); hold on
plot_rawData_operantMatching(fig_rawData, sessionData, unCorrectedBlockSwitch, tMax);

% P(actionChosen|model) plot
subplot(fig_probAction); hold on

probChosenChoice_qLearning = zeros(length(probChoice_qLearning), 1);
probChosenChoice_directActor = zeros(length(probChoice_directActor), 1);
probChosenChoice_acf = zeros(length(probChoice_acf), 1);
probChosenChoice_dynamicThompson = zeros(length(probChoice_dynamicThompson), 1);
probChosenChoice_stackedProb = zeros(length(probChoice_stackedProb), 1);
probChosenChoice_ph = zeros(length(probChoice_ph), 1);
for i = 1:length(probChoice_qLearning)
    tmp = probChoice_qLearning(i, :);
    probChosenChoice_qLearning(i) = tmp(logical(choice(i, :)));
    
    tmp = probChoice_directActor(i, :);
    probChosenChoice_directActor(i) = tmp(logical(choice(i, :)));
    
    tmp = probChoice_acf(i, :);
    probChosenChoice_acf(i) = tmp(logical(choice(i, :)));
    
    tmp = probChoice_dynamicThompson(i, :);
    probChosenChoice_dynamicThompson(i) = tmp(logical(choice(i, :)));
    
    tmp = probChoice_stackedProb(i, :);
    probChosenChoice_stackedProb(i) = tmp(logical(choice(i, :)));
    
    tmp = probChoice_ph(i, :);
    probChosenChoice_ph(i) = tmp(logical(choice(i, :)));
end
plot(xvals, probChosenChoice_qLearning,'Color', qLearning_color,'linewidth',1)
plot(xvals, probChosenChoice_directActor, 'Color', directActor_color,'linewidth',1)
plot(xvals, probChosenChoice_acf, 'Color', acf_color, 'linewidth', 1);
plot(xvals, probChosenChoice_dynamicThompson, 'Color', dynamicThompson_color,'linewidth',1)
plot(xvals, probChosenChoice_stackedProb, 'Color', stackedProb_color,'linewidth',1)
plot(xvals, probChosenChoice_ph, 'Color', ph_color, 'linewidth', 1)
plot([1 tMax], [0.5 0.5],'k--')
xlabel('Trials')
ylabel('Probability of Action Taken')
xlim([1 tMax])
ylim([0 1])   

% legend('Q-Learning', 'Direct Actor', 'Actor Critic', 'Dynamic Thompson', 'Stacked Probability', 'PRO')

% P(actionToChoose|model) plot
subplot(fig_compareModelAndSmoothBehavior); hold on

normKern = normpdf(-15:15,0,2);
normKern = normKern / sum(normKern);

smoothChoice = conv(behavStruct.allChoices, normKern, 'same');
smoothReward = conv(behavStruct.allRewards, normKern, 'same');
plot(xvals, smoothChoice, 'k', 'linewidth', 2);
plot(xvals, smoothReward, '--', 'Color', [100 100 100]./255, 'linewidth', 2)
xlabel('Trials')
ylabel('<-- Left       Right -->')
xlim([1 tMax])
ylim([-1 1])

plot(xvals, probChoice_qLearning(:, 1)*2 - 1, 'Color', qLearning_color, 'linewidth', 1)
plot(xvals, probChoice_directActor(:, 1)*2 - 1, 'Color', directActor_color, 'linewidth', 1)
plot(xvals, probChoice_acf(:, 1)*2 - 1, 'Color', acf_color, 'linewidth', 1);
plot(xvals, probChoice_dynamicThompson(:, 1)*2 - 1, 'Color', dynamicThompson_color, 'linewidth', 1)
plot(xvals, probChoice_stackedProb(:, 1)*2 - 1, 'Color', stackedProb_color, 'linewidth', 1)
plot(xvals, probChoice_ph(:, 1)*2 - 1, 'Color', ph_color, 'linewidth', 1)

% legend('Choices','Rewards', 'Q-Learning', 'Direct Actor', 'Actor Critic Fukai', 'Dynamic Thompson','Stacked Probability', 'PRO')


% Model selection plots
modelLabels = {};
for i = 1:length(qLearningModels)
    modelLabels = [modelLabels {['q_' int2str(i)]}];
end
for i = 1:length(directActorModels)
    modelLabels = [modelLabels {['da_' int2str(i)]}];
end
for i = 1:length(directActorWithValueModels)
    modelLabels = [modelLabels {['dav_' int2str(i)]}];
end
for i = 1:length(acfModels)
    modelLabels = [modelLabels {['acf_' int2str(i)]}];
end
for i = 1:length(dynamicThompsonModels)
    modelLabels = [modelLabels {['dts_' int2str(i)]}];
end
for i = 1:length(stackedProbModels)
    modelLabels = [modelLabels {['sp_' int2str(i)]}];
end
for i = 1:length(phModels)
    modelLabels = [modelLabels {['ph_' int2str(i)]}];
end

separation = [length(qLearningModels) length(directActorModels) length(directActorWithValueModels) ...
              length(acfModels) length(dynamicThompsonModels) length(stackedProbModels) length(phModels)];
separation = cumsum(separation) + 0.5; % lines to divide models


width = 0.8;
subplot(fig_LH); hold on
allLH = [qLH daLH davLH acfLH dtsLH spLH phLH];
b(1) = bar(allLH, width);
LHrange = max(allLH) - min(allLH);
ylim([min(allLH) - 0.3*LHrange, max(allLH) + 0.3*LHrange])
ylabel('Log Likelihood Score')
plot([separation; separation],[ones(1, length(separation))*-10000; ones(1, length(separation))*10000],'k--')

subplot(fig_BICscores); hold on
allBIC = [qBIC daBIC davBIC acfBIC dtsBIC spBIC phBIC];
b(2) = bar(allBIC, width);
BICrange = max(allBIC) - min(allBIC);
ylim([min(allBIC) - 0.3*BICrange, max(allBIC) + 0.3*BICrange])
ylabel('BIC Score')
plot([separation; separation],[ones(1, length(separation))*-10000; ones(1, length(separation))*10000],'k--')

set(fig_LH,'XTickLabel',modelLabels, 'XTick',1:length(allLH), 'TickLabelInterpreter', 'none')
set(fig_BICscores,'XTickLabel',modelLabels, 'XTick',1:length(allBIC), 'TickLabelInterpreter', 'none')

% plot best model for each class in red
all_bestInd = [best_qLearningModel best_directActorModel best_directActorWithValueModel best_acfModel best_dynamicThompsonModel best_stackedProbModel best_phModel];
modelLen = [length(qBIC) length(daBIC) length(davBIC) length(acfBIC) length(dtsBIC) length(spBIC) length(phBIC)];
modelLen = cumsum(modelLen);
all_bestInd(2:end) = all_bestInd(2:end) + modelLen(1:end-1);
subplot(fig_LH)
bar(all_bestInd, allLH(all_bestInd), width/min(diff(sort(all_bestInd))), 'r')
subplot(fig_BICscores)
bar(all_bestInd, allBIC(all_bestInd), width/min(diff(sort(all_bestInd))), 'r')


% plot overall best model in green
[~, bestMod_Ind] = min(allBIC);
subplot(fig_LH)
bar(bestMod_Ind, allLH(bestMod_Ind), width, 'g')
subplot(fig_BICscores)
bar(bestMod_Ind, allBIC(bestMod_Ind), width, 'g')


allSubplots = [fig_compareModelAndSmoothBehavior; fig_rawData; fig_probAction];
set(allSubplots,'Ticklength', [0 0])

for j = 1:length(allSubplots)
    subplot(allSubplots(j))
    for i = 1:length(unCorrectedBlockSwitch)
        bs_loc = unCorrectedBlockSwitch(i);
        plot([bs_loc bs_loc],[-1 1],'--','Color',[30 144 255]./255);
    end
end
subplot(fig_probAction)
legend('Q-Learning', 'Direct Actor', 'Actor Critic', 'Dynamic Thompson', 'Stacked Probability', 'Pearce-Hall')
subplot(fig_compareModelAndSmoothBehavior)
legend('Choices','Rewards', 'Q-Learning', 'Direct Actor', 'Actor Critic Fukai', 'Dynamic Thompson','Stacked Probability', 'Pearce-Hall')

% save figure
set(gcf, 'Position', get(0,'Screensize'))

if saveFigFlag == 1
    savepath = [root out.animalName sep out.sessionFolder sep  'figures' sep];
    if isempty(dir(savepath))
        mkdir(savepath)
    end

    % sessionFolder(~(sessionFolder==sep)) removes the separator (/ or \) in the filename before saving
    saveFigurePDF(gcf,[savepath sep out.saveFigName '_modelFits'])
    
    savePathModels = [root out.animalName sep out.sessionFolder sep 'sorted' sep out.saveFigName '_modelFits.mat'];
    save(savePathModels, 'models');
end