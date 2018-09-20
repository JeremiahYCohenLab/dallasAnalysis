function [LH, probChoice, Q, pe] = qLearningModel_2params(startValues, choice, outcome)

alphaLearn = startValues(1);
beta = startValues(2);

trials = length(choice);
Q = zeros(trials,2);

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice       
        Q(t+1, 2) = Q(t, 2);
        Q(t+1, 1) = Q(t, 1) + alphaLearn * (outcome(t, 1) - Q(t, 1));
        pe(t) = outcome(t, 1) - Q(t, 1);
    else % left choice
        Q(t+1, 1) = Q(t, 1);
        Q(t+1, 2) = Q(t, 2) + alphaLearn * (outcome(t, 2) - Q(t, 2));
        pe(t) = outcome(t, 2) - Q(t, 2);
    end
end

% Call softmax  rule

probChoice = logistic([beta*(Q(:, 1)-Q(:, 2)), ...
                       beta*(Q(:, 2)-Q(:, 1))]);

% To calculate likelihood:
LH = likelihood(choice,probChoice);
end