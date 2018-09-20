function [optimalParams, LH, exitFl, CIbands, probChoice, actions] = directActor(fileOrFolder)


% Pull out reward and choice data
if ischar(fileOrFolder)
    fileOrFolder = {fileOrFolder};
end
outcome = [];
choice = [];
sessionChangeInd = [];
for i = 1:length(fileOrFolder)    
    [behSessionData, blockSwitch] = loadBehavioralData(fileOrFolder{i});
    behavStruct = parseBehavioralData(behSessionData, blockSwitch);

    outcome = [outcome; abs([behavStruct.allReward_R; behavStruct.allReward_L])'];
    choice = [choice; abs([behavStruct.allChoice_R; behavStruct.allChoice_L])'];
    sessionChangeInd = [sessionChangeInd length(outcome)];
end

% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

alphaAction_range = [0 1];
alphaForgetAction_range = [0 1];
beta_range = [0 double(intmax)];
bias_range = [-5 5];


startValues = [0.5, 0.5, 10, 0];
A=[eye(length(startValues)); -eye(length(startValues))];
b=[ alphaAction_range(2);  alphaForgetAction_range(2);  beta_range(2);  bias_range(2);
   -alphaAction_range(1); -alphaForgetAction_range(1); -beta_range(1); -bias_range(1)];


[optimalParams(1,:), LH(1,:), exitFl(1,:), ~, ~, ~, hess] = fmincon(@directActorLearning,startValues,...
            A,b,[],[],[],[],[],options,choice,outcome,sessionChangeInd);
        
CIbands = sqrt(diag(inv(hess)))'*1.96;
        

[probChoice, actions] = plotFigure(fileOrFolder, optimalParams,choice,outcome,sessionChangeInd,behavStruct,blockSwitch);


function [LH, P, actions] = directActorLearning(startValues, choice, outcome, sessionChangeInd)

    
alphaAction = startValues(1);
alphaForgetAction = startValues(2);
beta = startValues(3);
bias = startValues(4);

trials = length(choice);
actions = zeros(trials,2);
P = zeros(trials,2);

% Call learning rule
t = 1;
P(t,:) = logistic([beta*(actions(t,1) - actions(t,2)) + bias, ...
                   beta*(actions(t,2) - actions(t,1)) - bias]);
for t = 1 : (trials-1)
    if any(t == sessionChangeInd) % if the session has changed
        actions(t,:) = 0; % restart action parameters
    end
    
    if choice(t, 1) == 1 % right choice
        
        actions(t+1, 1) = alphaForgetAction*actions(t, 1) + alphaAction*(1 - P(t,1))*outcome(t,1);
        actions(t+1, 2) = alphaForgetAction*actions(t, 2) - alphaAction*P(t,2)*outcome(t,1);
    else % left choice
       
        actions(t+1, 2) = alphaForgetAction*actions(t, 2) + alphaAction*(1 - P(t,2))*outcome(t,2);
        actions(t+1, 1) = alphaForgetAction*actions(t, 1) - alphaAction*P(t,1)*outcome(t,2);
    end
    P(t+1,:) = logistic([beta*(actions(t+1,1) - actions(t+1,2)) + bias, ...
                         beta*(actions(t+1,2) - actions(t+1,1)) - bias]);
end

% To calculate likelihood:
Pf = choice.* P;
nonzeros = find(Pf ~= 0);
Pf = Pf(nonzeros);
% LH = -1 * prod(Pf,1);
LH = -1 * sum(log(Pf));
end


function [probChoice, actions] = plotFigure(fileOrFolder, optimalParams, choice, outcome, sessionChangeInd, behavStruct, blockSwitch)
    
    [~,probChoice, actions] = directActorLearning(optimalParams, choice, outcome, sessionChangeInd);

    normKern = normpdf(-15:15,0,2);
    normKern = normKern / sum(normKern);
    xVals = (1:(length(normKern) + length(behavStruct.allChoices) - 1)) - round(length(normKern)/2);

    figure; 
    subplot(2,1,1); hold on
    plot(xVals, conv(behavStruct.allChoices,normKern),'k','linewidth',2);
    plot(xVals, conv(behavStruct.allRewards,normKern),'--','Color',[100 100 100]./255,'linewidth',2)
    xlabel('Trials')
    ylabel('<-- Left       Right -->')
    xlim([1 length(behavStruct.allChoice_R)])
    ylim([-1 1])

    plot(probChoice(:,1)*2-1,'-','Color',[148,0,211]./255,'linewidth',1)

    legend('Choices','Rewards', ...
            sprintf('%s\n%s: %s\n%s: %s\n%s: %s\n%s: %s\n%s: %s', 'Q-learning Model',....
            'alphaAction',num2str(optimalParams(1),2), ...
            'alphaForgetAction',num2str(optimalParams(2),2), ...
            'beta',num2str(optimalParams(3),2), ...
            'bias',num2str(optimalParams(4),2)))

    for i = 1:length(blockSwitch)
        bs_loc = blockSwitch(i);
        plot([bs_loc bs_loc],[-1 1],'--','Color',[30 144 255]./255)
    end

    title(fileOrFolder)

    subplot(2,1,2); hold on
    plot(actions,'linewidth',2); 
    ylabel('Actions')
    xlabel('Trials')
    legend('Right','Left')
    xlim([1 length(behavStruct.allChoice_R)])
    ylim([min(min(actions)) max(max(actions))])

    for i = 1:length(blockSwitch)
        bs_loc = blockSwitch(i);
        plot([bs_loc bs_loc],[0 1],'--','Color',[30 144 255]./255)
    end
end

end