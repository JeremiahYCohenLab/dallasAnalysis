% function analyzeOperantMatching_20160722
sessionToImport = '/Volumes/bbari1/BB039/mBB039d20160819/mBB039d20160819.asc';

tStart = 1;

sessionText = importMatchingData(sessionToImport);

stringList = {'Trial ', ...
              'No lick period', ...
              'CS PLUS', ...
              'CS MINUS', ...
              'Odor on', ...
              'Odor off', ...
              'ITI - ', ...
              'WATER L DELIVERED', ...
              'WATER L NOT DELIVERED', ...
              'WATER R DELIVERED', ...
              'WATER R NOT DELIVERED', ...
              'L:', ...
              'R:'};

 outputs = [];
for currStr = 3:length(stringList)
    tempInd = strfind(sessionText, stringList{currStr});
    strInds = find(not(cellfun('isempty', tempInd)));
    tempExtract = regexp(sessionText(strInds), ': ', 'split');
    allStr = [tempExtract{:}];
    outputs{currStr} = cellfun(@str2num,allStr(2:2:end));
end

% R = 1, L = -1
allLicks = sort([outputs{8} outputs{9} outputs{10} outputs{11}]);
allLickDecisions = allLicks;
allLickDecisions(ismember(allLickDecisions, outputs{10})) = 1;  % all R licks = 1
allLickDecisions(ismember(allLickDecisions, outputs{11})) = 1;
allLickDecisions(ismember(allLickDecisions, outputs{8})) = -1;  % all L licks = -1
allLickDecisions(ismember(allLickDecisions, outputs{9})) = -1;

allRewards = allLicks;
allRewards(ismember(allRewards, outputs{10})) = 1; % all R rewards = 1
allRewards(ismember(allRewards, outputs{8})) = -1; % all L rewards = 0
allRewards(allRewards ~=1 & allRewards ~=-1) = 0;

figure; subplot(3,1,1); hold on; 
plot(smooth(allLickDecisions(tStart:end),10),'k','linewidth',2)
plot(smooth(allRewards(tStart:end),10),'b','linewidth',2)
ylim([-1 1])

rwdMatx = [];
for i = 1:10
    rwdMatx(i,:) = [NaN(1,i) allRewards(tStart:end-i)];
end

lickMatx = [];
for i = 1:10
    lickMatx(i,:) = [NaN(1,i) allLickDecisions(tStart:end-i)];
end
lm = fitlm([rwdMatx]', allLickDecisions(tStart:end),'linear');

subplot(3,1,2); plot(lm.Coefficients.Estimate(2:end))
ylim([-0.1 0.8])

choices = [];
rewards = [];
chunkSize = 30;
for choiceChunk = 1:chunkSize:length(allLickDecisions)-chunkSize
    choices(((chunkSize-1)+choiceChunk)/chunkSize) = sum(allLickDecisions(choiceChunk:choiceChunk+chunkSize-1) == 1)/chunkSize;
    rewards(((chunkSize-1)+choiceChunk)/chunkSize) = sum(allRewards(choiceChunk:choiceChunk+chunkSize-1) == 1)/ ...
        sum(allRewards(choiceChunk:choiceChunk+chunkSize-1) == 1 | allRewards(choiceChunk:choiceChunk+chunkSize-1) == -1);
end

subplot(3,1,3); plot(rewards, choices, 'o'); xlabel('Rewards (R > L)'); ylabel('Choices (R > L)');
xlim([0 1]); ylim([0 1])

%%
changeChoice = [false abs(diff(allLickDecisions)) > 0];
changeHistogram = [];
for i = find(changeChoice == 1)
    if allLickDecisions(i) == 1 % if a right lick
        temp = 0;
        goBack = 1;
        while (i - goBack > 0) && allLickDecisions(i-goBack) == -1 && allRewards(i-goBack) == 0 % if previous trial was a L lick AND resulted in no reward
            temp = temp + 1;
            goBack = goBack + 1;
        end
        changeHistogram = [changeHistogram temp];
    elseif allLickDecisions(i) == -1 %if a left lick
        temp = 0;
        goBack = 1;
        while (i - goBack > 0) && allLickDecisions(i-goBack) == 1 && allRewards(i-goBack) == 0 % if previous trial was a R lick AND resulted in no reward
            temp = temp + 1;
            goBack = goBack + 1;
        end
        changeHistogram = [changeHistogram temp];
    end
end

figure; histogram(changeHistogram,0:1:50,'Normalization','probability')

%%
changeChoice = [false abs(diff(allLickDecisions)) > 0];
changeHistogram_LtoR = [];
changeHistogram_RtoL = [];
for i = find(changeChoice == 1)
    if allLickDecisions(i) == 1 % if a right lick
        temp = 0;
        goBack = 1;
        while (i - goBack > 0) && allLickDecisions(i-goBack) == -1 && allRewards(i-goBack) == 0 % if previous trial was a L lick AND resulted in no reward
            temp = temp + 1;
            goBack = goBack + 1;
        end
        changeHistogram_LtoR = [changeHistogram_LtoR temp];
    elseif allLickDecisions(i) == -1 %if a left lick
        temp = 0;
        goBack = 1;
        while (i - goBack > 0) && allLickDecisions(i-goBack) == 1 && allRewards(i-goBack) == 0 % if previous trial was a R lick AND resulted in no reward
            temp = temp + 1;
            goBack = goBack + 1;
        end
        changeHistogram_RtoL = [changeHistogram_RtoL temp];
    end
end

figure; hold on;
histogram(changeHistogram_LtoR,0:50,'Normalization','probability')
histogram(changeHistogram_RtoL,0:50,'Normalization','probability')
legend('L -> R','R -> L')

%% generate figures
figure; hold on;
bS = [1    38    49   116   130   169   181   209];
avgRwdSlope = [0.4493    2.7350    0.3326    5.3509    0.5795    4.7876    0.7161]
plot(cumsum(allLickDecisions==-1), cumsum(allLickDecisions==1),'linewidth',2,'Color',[30,144,255]/255);

plot([bS(1) bS(2)],[bS(1) bS(2)]*avgRwdSlope(1),'k','linewidth',2);
plot([bS(2) bS(3)],[bS(2) bS(3)]*avgRwdSlope(2)-80,'k','linewidth',2);
plot([bS(3) bS(4)],[bS(3) bS(4)]*avgRwdSlope(3)+50,'k','linewidth',2);
plot([bS(4) bS(5)],[bS(4) bS(5)]*avgRwdSlope(4)-530,'k','linewidth',2);
plot([bS(5) bS(6)],[bS(5) bS(6)]*avgRwdSlope(5)+75,'k','linewidth',2);
plot([bS(6) bS(7)],[bS(6) bS(7)]*avgRwdSlope(6)-630,'k','linewidth',2);
plot([bS(7) bS(8)],[bS(7) bS(8)]*avgRwdSlope(7)+135,'k','linewidth',2)
legend('Choice Ratio','Income Ratio','location','best')
xlabel('Cumulative Left Choices'); ylabel('Cumulative Right Choices')
% xlim([0 300]); ylim([0 300])

%%
figure; hold on;
plot(conv(diff((cumsum(allLickDecisions == -1) - cumsum(allLickDecisions == 1))),normKern),'linewidth',2,'Color',[30,144,255]/255)
plot(conv(diff((cumsum(allRewards == -1) - cumsum(allRewards == 1))),normKern),'k','linewidth',2)

%%
figure; hold on;
bS = [1 55 110 208 277 337 435 509];
normKern = normpdf(-15:15,0,6);
normKern = normKern / sum(normKern);
normKern = normKern(15:end);

lickSlope = 57.2958*atan(diff(conv(cumsum(allLickDecisions == 1),normKern))./diff(conv(cumsum(allLickDecisions == -1),normKern)));
rwdSlope = 57.2958*atan(diff(conv(cumsum(allRewards == 1),normKern))./diff(conv(cumsum(allRewards == -1),normKern)));
plot(lickSlope,'linewidth',1.5,'Color',[30,144,255]/255)
plot(rwdSlope,'k','linewidth',1.5)
avgRwdSlope = [];
% for i = 2:length(bS)
%     plot([bS(i-1) bS(i)], [mean(lickSlope(bS(i-1):bS(i))) mean(lickSlope(bS(i-1):bS(i)))],'linewidth',4,'Color',[30,144,255]/255);
%     plot([bS(i-1) bS(i)], [mean(rwdSlope( bS(i-1):bS(i))) mean(rwdSlope( bS(i-1):bS(i)))],'k','linewidth',4);
%     avgRwdSlope = [avgRwdSlope tand(mean(rwdSlope( bS(i-1):bS(i))))];
% end

xlim([0 510])
ylim([0 90])
legend('Instantaneous Choice Ratio','Instantaneous Income Ratio','location','bestoutside')
xlabel('Choice number')
ylabel('Slope (degrees)')

%% reward history
a = 0.3689;
b = -0.2973;
x = 0:.1:50; 
expFit = a*exp(b*x);
figure; hold on; xlim([1 50])
plot(lm.Coefficients.Estimate(2:end)/sum(lm.Coefficients.Estimate(2:end)),'ko');
plot(x,expFit,'k','linewidth',2)
xlabel('Time Since Rewards (trials)'); ylabel('Relative Filter Weight');
legend('Mouse','Best Fit');