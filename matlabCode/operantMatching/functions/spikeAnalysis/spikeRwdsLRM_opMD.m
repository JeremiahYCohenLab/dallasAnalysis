function [mdlStruct, betaMaxTau, betaPreTau] = spikeRwdsLRM_opMD(cellFile, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('timeMax', 121000);
p.addParameter('timeBins', 12);
p.addParameter('tMax', 12);
p.addParameter('trialFlag', 0);
p.addParameter('plotFlag', 1);
p.parse(varargin{:});

[root, sep] = currComputer();
savePath = ['C:\Users\cooper_PC\Desktop\analysis\rwdsOnSpikesLRM\'];

[revForFlag, cellList, ~] = xlsread(cellFile, category);
cellList = cellList(2:end,:);

smoothWin = 250;
tb = 1.5;
tf = 5;
trialBeg = tb*1000;
CSoff = tb*1000 + 500;
time = -1000*tb:1000*tf;

if p.Results.trialFlag
    binArray = [1:p.Results.tMax];
    numBins = p.Results.tMax;
else
    timeBinSize = (p.Results.timeMax - 1000)/p.Results.timeBins;
    binArray = [(1:p.Results.timeBins)*timeBinSize/1000];
    numBins = p.Results.timeBins;
end
    
mdlStruct = struct;
betaMaxTau = [];
betaPreTau = [];
plotMaxFit = 0;
plotPreFit = 0;

for i = 1:length(cellList)
    i=i+0
    if i==1 | isempty(regexp(cellList{i,2}, cellList{i-1,2}))
        sessionName = cellList{i,2};
        [animalName] = strtok(sessionName, 'd');
        animalName = animalName(2:end);

        if isstrprop(sessionName(end), 'alpha')
            sortedFolderLocation = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
        else
            sortedFolderLocation = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
        end
        sortedFolder = dir(sortedFolderLocation);

    %     if intanFlag{i} == 1
    %         if any(~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'))) == 1
    %             sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'));
    %             load([sortedFolderLocation sortedFolder(sessionDataInd).name])
    %         else
    %             [sessionData] = generateSessionData_intan_operantMatching(sessionName);
    %         end
    %     else
            if any(~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'))) == 1
                sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'));
                load([sortedFolderLocation sortedFolder(sessionDataInd).name])
            else
                [sessionData] = generateSessionData_nL_operantMatching(sessionName);
            end
    %    end
        [s] = behAnalysisNoPlot_opMD(sessionName, 'revForFlag', revForFlag(i), 'timeMax', p.Results.timeMax,...
            'timeBins', p.Results.timeBins, 'tMax', p.Results.tMax);
    end
    

    %% Sort all spikes into a raster-able matrix
    
    cellName = cellList{i,1};
    spikeFields = fields(sessionData);
    cellInd(i) = find(~cellfun(@isempty,strfind(spikeFields,cellName)));

    allTrial_spike = {};
    for j = 1:length(sessionData)
        if j == 1
            prevTrial_spike = [];
            currTrial_lick = [];
        else
            prevTrial_spikeInd = sessionData(j-1).(spikeFields{cellInd(i)}) > sessionData(j-1).trialEnd-tb*1000;
            prevTrial_spike = sessionData(j-1).(spikeFields{cellInd(i)})(prevTrial_spikeInd) - sessionData(j).CSon;
        end

        currTrial_spikeInd = sessionData(j).(spikeFields{cellInd(i)}) < sessionData(j).CSon+tf*1000;
        currTrial_spike = sessionData(j).(spikeFields{cellInd(i)})(currTrial_spikeInd) - sessionData(j).CSon;

        allTrial_spike{i,j} = [prevTrial_spike currTrial_spike];
    end

    % sometimes no licks/spikes are considered 1x0 and sometimes they are []
    % plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
    allTrial_spike(cellfun(@isempty,allTrial_spike)) = {zeros(1,0)}; 

    for j = 1:length(s.behSessionData)
        trialDurDiff(j) = (s.behSessionData(j).trialEnd - s.behSessionData(j).CSon)- tf*1000;
    end

    allTrial_spikeMatx = NaN(length(sessionData),length(time));
    for j = 1:length(allTrial_spike)
        tempSpike = allTrial_spike{j};
        tempSpike = tempSpike + tb*1000; % add this to pad time for SDF
        allTrial_spikeMatx(j,tempSpike) = 1;
        if trialDurDiff(j) < 0
            allTrial_spikeMatx(j, isnan(allTrial_spikeMatx(j, 1:end+trialDurDiff(j)))) = 0;  %converts within trial duration NaNs to 0's
        else
            allTrial_spikeMatx(j, isnan(allTrial_spikeMatx(j,:))) = 0;
        end
        if sum(allTrial_spikeMatx(j,:)) == 0     %if there is no spike data for this trial, don't count it
            allTrial_spikeMatx(j,:) = NaN;
        end
    end
    
    %% plot effect of previous reward
    prevRwd_spikeMatx = [];
    prevNoRwd_spikeMatx = [];

    for j = 2:length(s.responseInds)
        if s.allRewardsBinary(j-1) == 1
            prevRwd_spikeMatx = [prevRwd_spikeMatx; allTrial_spikeMatx(s.responseInds(j),:)];
        else
            prevNoRwd_spikeMatx = [prevNoRwd_spikeMatx; allTrial_spikeMatx(s.responseInds(j),:)];
        end
    end

    mySDF_prevRwd = fastsmooth(nanmean(prevRwd_spikeMatx, 1)*1000, smoothWin)';
    mySDF_prevNoRwd = fastsmooth(nanmean(prevNoRwd_spikeMatx, 1)*1000, smoothWin)';
    mySD_prevRwd = fastsmooth(nanstd(prevRwd_spikeMatx, 1)*1000 / sqrt(size(prevRwd_spikeMatx,2)), smoothWin)';
    mySD_prevNoRwd = fastsmooth(nanstd(prevNoRwd_spikeMatx, 1)*1000 / sqrt(size(prevRwd_spikeMatx,2)), smoothWin)';
   
    
    if p.Results.plotFlag
        fig = figure; hold on;
        subplot(1,3,1)
        plot(mySDF_prevNoRwd(smoothWin:end-smoothWin),'-','Color', [0.5 0 0.8],'LineWidth',2); hold on
        plot(mySDF_prevRwd(smoothWin:end-smoothWin),'b','LineWidth',2)
        legend('prev no rwd', 'prev rwd');
    end
   
    %% find features of spike rate on each trial    
    
    preCSspikeCount = NaN(1,length(allTrial_spike));
    postCSspikeCount = NaN(1,length(allTrial_spike));
    maxFRcs = NaN(1,length(allTrial_spike));
    minFRcs = NaN(1,length(allTrial_spike));
    maxFRtime = NaN(1,length(allTrial_spike));
    
    for j = 1:length(allTrial_spike)
        if ~isempty(allTrial_spikeMatx(j))
            preCSspikeCount(j) = sum(allTrial_spikeMatx(j, 1:trialBeg));              %find total spikes before CS on
            postCSspikeCount(j) = sum(allTrial_spikeMatx(j, trialBeg:CSoff));

            spikeTemp = fastsmooth(allTrial_spikeMatx(j,:)*1000, smoothWin, 3);         %smooth raw spikes to find features of spike rate
            maxFRcs(j) = max(spikeTemp(trialBeg:CSoff));
            minFRcs(i,j) = min(spikeTemp(trialBeg:CSoff));
            if ~isnan(maxFRcs(j))
                maxFRtime(j) = find(spikeTemp == max(spikeTemp(trialBeg:CSoff)), 1);
            else
                maxFRtime(j) = NaN;
            end
        end
    end
    
    nameTmp = [sessionName '_' cellName];
    if p.Results.trialFlag
        mdlStruct.(nameTmp).maxFR = fitlm([s.rwdMatx]', maxFRcs(s.responseInds));
        mdlStruct.(nameTmp).preCS = fitlm([s.rwdMatx]', preCSspikeCount(s.responseInds));
    else
        mdlStruct.(nameTmp).maxFR = fitlm([s.rwdTimeMatx]', maxFRcs(s.responseInds));
        mdlStruct.(nameTmp).preCS = fitlm([s.rwdTimeMatx]', preCSspikeCount(s.responseInds));
    end
    
    if mdlStruct.(nameTmp).maxFR.Coefficients.pValue(2) < 0.05
        ft = fittype( 'a*exp((-1/b)*x)', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares'); opts.Display = 'Off'; opts.StartPoint = [-3 20];
        expFitMax = fit(binArray', mdlStruct.(nameTmp).maxFR.Coefficients.Estimate(2:numBins+1), ft, opts );
        fitCoeffsMax = coeffvalues(expFitMax); 
        betaMaxTau = [betaMaxTau; {nameTmp} fitCoeffsMax(2)];
        plotMaxFit = 1;
    end
    
    
    if mdlStruct.(nameTmp).preCS.Coefficients.pValue(2) < 0.05
        ft = fittype( 'a*exp((-1/b)*x)', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares'); opts.Display = 'Off'; opts.StartPoint = [-3 20];
        expFitPre = fit(binArray', mdlStruct.(nameTmp).preCS.Coefficients.Estimate(2:numBins+1), ft, opts );
        fitCoeffsPre = coeffvalues(expFitPre); 
        betaPreTau = [betaPreTau; {nameTmp} fitCoeffsPre(2)];
        plotPreFit = 1;
    end
        
    if p.Results.plotFlag
        subplot(1,3,2); hold on;
        relevInds = 2:p.Results.timeBins+1;
        coefVals = mdlStruct.(nameTmp).maxFR.Coefficients.Estimate(relevInds);
        CIbands = coefCI(mdlStruct.(nameTmp).maxFR);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar(binArray,coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
        ylabel('\beta Coefficient')
        if p.Results.trialFlag
            xlabel('reward n trials back')
            title('rewards in trials on max fr')
        else
            xlabel('reward n seconds back')
            title('rewards in time on max fr')
        end
        xlim([0 (binArray(end) + 1)])
        if plotMaxFit
            y = fitCoeffsMax(1)*exp(-(1/fitCoeffsMax(2))*binArray);
            plot(binArray, y, '-c', 'LineWidth', 2)
            plotMaxFit = 0;
            xl = get(gca, 'XLim'); yl = get(gca, 'YLim');
            text('Position', [xl(2) yl(2) 0], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top',...
                'string', ['tau: ' num2str(fitCoeffsMax(2))])
        end

        subplot(1,3,3); hold on;
        coefVals = mdlStruct.(nameTmp).preCS.Coefficients.Estimate(relevInds);
        CIbands = coefCI(mdlStruct.(nameTmp).preCS);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar(binArray,coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
        ylabel('\beta Coefficient')
        if p.Results.trialFlag
            xlabel('reward n trials back')
            title('rewards in trials on max fr')
        else
            xlabel('reward n seconds back')
            title('rewards in time on max fr')
        end
        xlim([0 (binArray(end) + 1)])
        cellNameTmp = erase(cellName, '_');
        suptitle([sessionName ': ' cellNameTmp])
        
        if plotPreFit
            y = fitCoeffsPre(1)*exp(-(1/fitCoeffsPre(2))*binArray);
            plot(binArray, y, '-c', 'LineWidth', 2)
            plotPreFit = 0;
            xl = get(gca, 'XLim'); yl = get(gca, 'YLim');
            text('Position', [xl(2) yl(2) 0], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top',...
                'string', ['tau: ' num2str(fitCoeffsPre(2))])
        end
    end
    
    set(fig, 'Position', get(0,'Screensize'))
    if p.Results.trialFlag
        saveFigurePDF(fig, [savePath 'trial\' category sep sessionName '_' cellName]);
    else
        saveFigurePDF(fig, [savePath 'time\' category sep sessionName '_' cellName]);
    end
    close(fig);
    
end



    