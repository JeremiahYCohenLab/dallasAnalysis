function [LH, probChoice, Q, pe, rBar] = qLearningModel_5params_tForget_opponency(startValues, choice, outcome, ITI)

alphaNPE = startValues(1);
alphaPPE = startValues(2);
tForget = startValues(3);
beta = startValues(4);
v = startValues(5);

trials = length(choice);
Q = zeros(trials,2);
rBar = [0.4; zeros(trials-1,1)];

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        Q(t+1, 2) = Q(t, 2);
        rBar(t+1) = v*outcome(t, 1) + (1-v)*rBar(t);
        pe(t) = outcome(t, 1) - Q(t, 1)  - rBar(t);
        if pe(t) < 0
            Q(t+1, 1) = Q(t, 1) + alphaNPE * pe(t);
        else
            Q(t+1, 1) = Q(t, 1) + alphaPPE * pe(t);
        end
    else % left choice
        Q(t+1, 1) = Q(t, 1);
        rBar(t+1) = v*outcome(t, 2) + (1-v)*rBar(t);
        pe(t) = outcome(t, 2) - Q(t, 2) - rBar(t);
        if pe(t) < 0
            Q(t+1, 2) = Q(t, 2) + alphaNPE * pe(t);
        else
            Q(t+1, 2) = Q(t, 2) + alphaPPE * pe(t);
        end
    end
    Q(t+1, :) = Q(t+1, :)*exp(-tForget*ITI(t+1));
end

if choice(t, 1) == 1
    pe(trials) = outcome(end, 1) - Q(end, 1);
else
    pe(trials) = outcome(end, 2) - Q(end, 2);
end


% Call softmax  rule

probChoice = logistic([beta*(Q(:, 1)-Q(:, 2)), ...
                       beta*(Q(:, 2)-Q(:, 1))]);

% To calculate likelihood:
LH = likelihood(choice,probChoice);
end