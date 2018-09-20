function [LH, probChoice, Q, pe, R] = qLearningModel_4params_gradient(startValues, choice, outcome)

alphaNPE = startValues(1);
alphaPPE = startValues(2);
alphaR = startValues(3);
beta = startValues(4);

trials = length(choice);
Q = zeros(trials,2);
pe = zeros(trials,1);
R = [0.7; zeros(trials-1,1)];
probChoice = [0.5 0.5; nan(trials-1,2)];

% Call learning rule
for t = 1 : (trials-1)
    
    if choice(t, 1) == 1 % right choice
        pe(t) = outcome(t,1) - R(t);
        R(t+1) = R(t) + alphaR * pe(t);
        if pe(t) < 0
            Q(t+1, 1) = Q(t, 1) + alphaNPE * pe(t) * (1 - probChoice(t,1));
            Q(t+1, 2) = Q(t, 2) - alphaNPE * pe(t) * probChoice(t,2);
        else
            Q(t+1, 1) = Q(t, 1) + alphaPPE * pe(t) * (1 - probChoice(t,1));
            Q(t+1, 2) = Q(t, 2) - alphaPPE * pe(t) * probChoice(t,2);
        end
    else % left choice
        pe(t) = outcome(t,2) - R(t);
        R(t+1) = R(t) + alphaR * pe(t);
        if pe(t) < 0
            Q(t+1, 2) = Q(t, 2) + alphaNPE * pe(t) * (1 - probChoice(t,2));
            Q(t+1, 1) = Q(t, 1) - alphaNPE * pe(t) * probChoice(t,1);
        else
            Q(t+1, 2) = Q(t, 2) + alphaPPE * pe(t) * (1 - probChoice(t,2));
            Q(t+1, 1) = Q(t, 1) - alphaPPE * pe(t) * probChoice(t,1);
        end
    end
    probChoice(t+1, :) = logistic([beta*(Q(t+1, 1)-Q(t+1, 2)) beta*(Q(t+1, 2)-Q(t+1, 1))]);
end

if choice(t, 1) == 1
    pe(trials) = outcome(end, 1) - R(end);
else
    pe(trials) = outcome(end, 2) - R(end);
end

% To calculate likelihood:
LH = likelihood(choice,probChoice);
end