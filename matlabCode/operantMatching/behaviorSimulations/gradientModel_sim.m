function [Q, R] = gradientModel_sim(varargin)
%
%
% Simulate dynamic foraging task with temporally-forgetting Q learning model
% 
%
%task and model parameters
a = inputParser;
% default parameters if none given
a.addParameter('maxTrials', 1e3)
a.addParameter('blockLength', [20 35]);
a.addParameter('rwdProb', [90 10]);
a.addParameter('ITIparam', 0.3);
a.addParameter('bestParams', [0.1, 0.1, 5]);
a.parse(varargin{:});

alpha = a.Results.bestParams(1);
alphaR = a.Results.bestParams(2);
beta = a.Results.bestParams(3);


%initialize task class
p = RestlessBandit('RandomSeed',25,'blockLength', a.Results.blockLength,'maxTrials', a.Results.maxTrials, 'RewardProbabilities', a.Results.rwdProb);

% [left, right]; these are Q values going INTO that trial, before making a decision
Q = [0 0; NaN(a.Results.maxTrials-1, 2)];      % initialize Q values as 0
R = [0; NaN(a.Results.maxTrials-1, 1)];

% plot parameters
figure; 
rawData_plot = subplot(2, 1, 1); hold on; ylabel('<--- L          R --->');
qValue_plot = subplot(2, 1, 2); hold on; title('Q values'); xlabel('Trials'); ylabel('Q values');

for currT = 1:p.MaxTrials - 1
   
   % Select action
   pLeft = 1/(1 + exp(-beta*diff(Q(currT, :))));
   subplot(qValue_plot)
   title(sprintf('Probability of right choice %d%%', round(pLeft*100)))
   if binornd(1, pLeft) == 0 % left choice selected probabilistically
       p = p.inputChoice([1 0]);
       rpe = p.AllRewards(currT, 1) - R(currT);
       R(currT + 1) = R(currT) + alphaR * rpe;
       Q(currT + 1, 1) = Q(currT, 1) + alpha * rpe * (1 - pLeft);
       Q(currT + 1, 2) = Q(currT, 2) - alpha * rpe * (1 - pLeft);
   else
       p = p.inputChoice([0 1]);
       rpe = p.AllRewards(currT, 2) - R(currT);
       R(currT + 1) = R(currT) + alphaR * rpe;
       Q(currT + 1, 2) = Q(currT, 2) + alpha * rpe * (pLeft);
       Q(currT + 1, 1) = Q(currT, 1) - alpha * rpe * (pLeft);
   end
   

   subplot(rawData_plot); xlim([0 currT+1]); ylim([-1 1])
   if p.AllChoices(currT, 1) == 1 % left choice
       if p.AllRewards(currT, 1) == 1 % reward
           plot([currT, currT], [0 -1], 'k')
       else
           plot([currT, currT], [0 -0.5], 'k')
       end
   elseif p.AllChoices(currT, 2) == 1 % right choice
       if p.AllRewards(currT, 2) == 1 % reward
           plot([currT, currT], [0 1], 'k')
       else
           plot([currT, currT], [0 0.5], 'k')
       end
   end    
   
   subplot(qValue_plot); xlim([0 currT + 1]); ylim([-0.1 1.1]); ylim([0 max(max(Q))]);
   plot(Q(1:currT + 1, 1),'c')
   plot(Q(1:currT + 1, 2),'m')
   plot(R(1:currT + 1), '-', 'Color', [0.7 0 1]);
   legend('Left','Right')
end
suptitle('gradient algorithm simulated behavior');
end