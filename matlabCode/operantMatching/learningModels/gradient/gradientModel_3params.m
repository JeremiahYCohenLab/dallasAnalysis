function [LH, probChoice, Q, pe, R] = qLearningModel_3params_gradient(startValues, choice, outcome)

alpha = startValues(1);
alphaR = startValues(2);
beta = startValues(3);

trials = length(choice);
Q = zeros(trials,2);
pe = zeros(trials,1);
R = zeros(trials,1);
probChoice = [0.5 0.5; nan(trials-1,2)];

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        pe(t) = outcome(t, 1) - R(t);
        R(t+1) = R(t) + alphaR * pe(t);
        Q(t+1, 1) = Q(t, 1) + alpha * pe(t) * (1 - probChoice(t,1));
        Q(t+1, 2) = Q(t, 2) - alpha * pe(t) * probChoice(t,2);
    else % left choice
        pe(t) = outcome(t,2) - R(t);
        R(t+1) = R(t) + alphaR * pe(t);
        Q(t+1, 2) = Q(t, 2) + alpha * pe(t) * (1 - probChoice(t,2));
        Q(t+1, 1) = Q(t, 1) - alpha * pe(t) * probChoice(t,1);
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