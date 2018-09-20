function [LH, probChoice, Q_s, Q_l, pe_s, pe_l] = qLearningModel_dual(startValues, choice, outcome, ITI)

alphaNPE_s = startValues(1);
alphaPPE_s = startValues(2);
alphaForget_s = startValues(3);
beta = startValues(4);
ITIthresh = startValues(5);

alphaNPE_l = startValues(6);
alphaPPE_l = startValues(7);
alphaForget_l = startValues(8);

trials = length(choice);
Q_s = zeros(trials,2);
Q_l = zeros(trials,2);

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        Q_s(t+1, 2) = alphaForget_s*Q_s(t, 2);
        Q_l(t+1, 2) = alphaForget_l*Q_l(t, 2);
        pe_s(t) = outcome(t, 1) - Q_s(t, 1);
        pe_l(t) = outcome(t, 1) - Q_l(t, 1);
        if pe_s(t) < 0
            Q_s(t+1, 1) = alphaForget_s*Q_s(t, 1) + alphaNPE_s * pe_s(t);
        else
            Q_s(t+1, 1) = alphaForget_s*Q_s(t, 1) + alphaPPE_s * pe_s(t);
        end
        if pe_l(t) < 0
            Q_l(t+1, 1) = alphaForget_l*Q_l(t, 1) + alphaNPE_l * pe_l(t);
        else
            Q_l(t+1, 1) = alphaForget_l*Q_l(t, 1) + alphaPPE_l * pe_l(t);
        end
    else % left choice
        Q_s(t+1, 1) = alphaForget_s*Q_s(t, 1);
        Q_l(t+1, 1) = alphaForget_l*Q_l(t, 1);
        pe_s(t) = outcome(t, 2) - Q_s(t, 2);
        pe_l(t) = outcome(t, 2) - Q_l(t, 2);
        if pe_s(t) < 0
            Q_s(t+1, 2) = alphaForget_s*Q_s(t, 2) + alphaNPE_s * pe_s(t);
        else
            Q_s(t+1, 2) = alphaForget_s*Q_s(t, 2) + alphaPPE_s * pe_s(t);
        end
        if pe_l(t) < 0
            Q_l(t+1, 2) = alphaForget_l*Q_l(t, 2) + alphaNPE_l * pe_l(t);
        else
            Q_l(t+1, 2) = alphaForget_l*Q_l(t, 2) + alphaPPE_l * pe_l(t);
        end
    end 
end

if choice(t, 1) == 1
    pe_s(trials) = outcome(end, 1) - Q_s(end, 1);
    pe_l(trials) = outcome(end, 1) - Q_l(end, 1);
else
    pe_s(trials) = outcome(end, 2) - Q_s(end, 2);
    pe_l(trials) = outcome(end, 2) - Q_l(end, 2);
end

% Call softmax  rule

probChoice_s = logistic([beta*(Q_s(:, 1)-Q_s(:, 2)), ...
                       beta*(Q_s(:, 2)-Q_s(:, 1))]);
probChoice_l = logistic([beta*(Q_l(:, 1)-Q_l(:, 2)), ...
                       beta*(Q_l(:, 2)-Q_l(:, 1))]);

for t = 1:length(probChoice_s)
    if ITI(t) > ITIthresh
        probChoice(t,:) = probChoice_l(t,:);
    else
        probChoice(t,:) = probChoice_s(t,:);
    end
end

% To calculate likelihood:
LH = likelihood(choice,probChoice);
end