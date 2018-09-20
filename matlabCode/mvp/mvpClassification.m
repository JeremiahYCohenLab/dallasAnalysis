function [normFeatTbl, idx, D] = mvpClassification(plotFlag, labelFlag, numClusts, mvp)

if nargin < 1
    plotFlag = 1;   %flag for plotting results
end
if nargin < 2
    labelFlag = 0;    %flag for labeling data points on plot
end
if nargin < 3
    numClusts = 2;    %number of clusters specified for k-means clustering
end
if nargin < 4
    magno = load('C:\Users\cooper_PC\Desktop\mvp\clustering\magno.mat');        %load data
    parvo = load('C:\Users\cooper_PC\Desktop\mvp\clustering\parvo.mat');
    nac = load('C:\Users\cooper_PC\Desktop\mvp\clustering\nac.mat');
end

%get indeces for data groups
magnoInds = [1:length(magno.traceTbl.firstSpikeLat(:,1))];
parvoInds = [1:length(parvo.traceTbl.firstSpikeLat(:,1))] + magnoInds(end);
nacInds = [1:length(nac.traceTbl.firstSpikeLat(:,1))] + parvoInds(end);

%min-max normalization
% featMatx(:,1) = [magno.traceTbl.firstSpikeLat(:,1); parvo.traceTbl.firstSpikeLat(:,1); nac.traceTbl.firstSpikeLat(:,1)];
% featMatx(:,1) = (featMatx(:,1) - min(featMatx(:,1))) / (max(featMatx(:,1)) - min(featMatx(:,1)));
% featMatx(:,2) = [magno.traceTbl.spikeFreq(:,1); parvo.traceTbl.spikeFreq(:,1); nac.traceTbl.spikeFreq(:,1)];
% featMatx(:,2) = (featMatx(:,2) - min(featMatx(:,2))) / (max(featMatx(:,2)) - min(featMatx(:,2)));
% featMatx(:,3) = [magno.traceTbl.firstSpikeHW(:,1); parvo.traceTbl.firstSpikeHW(:,1); nac.traceTbl.firstSpikeHW(:,1)];
% featMatx(:,3) = (featMatx(:,3) - min(featMatx(:,3))) / (max(featMatx(:,3)) - min(featMatx(:,3)));

%zscore standardization
featMatx(:,1) = zscore([magno.traceTbl.firstSpikeLat(:,1); parvo.traceTbl.firstSpikeLat(:,1); nac.traceTbl.firstSpikeLat(:,1)]);
featMatx(:,2) = zscore([magno.traceTbl.spikeFreq(:,1); parvo.traceTbl.spikeFreq(:,1); nac.traceTbl.spikeFreq(:,1)]);
featMatx(:,3) = zscore([magno.traceTbl.firstSpikeHW(:,1); parvo.traceTbl.firstSpikeHW(:,1); nac.traceTbl.firstSpikeHW(:,1)]);


%k-means clustering analysis
[idx,C,sumd,D] = kmeans(featMatx(1:parvoInds(end),:), numClusts);


if plotFlag
    %determine plane that distinguishes between 2 generated clusters
    if numClusts == 2
        halfPoint = [C(2,:) + C(1,:)] / 2;
        normVec = C(2,:) - C(1,:);
        d = -[sum(halfPoint .* normVec)];
        [x, y] = meshgrid(-4:0.01:4);
        z = -1/normVec(3)*(normVec(1)*x + normVec(2)*y + d);
    end
    
    h = figure; 
    %h.Renderer = 'Painters';
    hold on; grid on;
    scatter3(featMatx(parvoInds,2), featMatx(parvoInds,3), featMatx(parvoInds,1), 'filled', 'MarkerFaceColor', [0.5 0.5 0.5]);
    scatter3(C(2,2), C(2,3), C(2,1), [100], 'Marker', '+', 'MarkerEdgeColor', [0.5 0.5 0.5]);
    scatter3(featMatx(magnoInds,2), featMatx(magnoInds,3), featMatx(magnoInds,1), 'filled', 'MarkerFaceColor', [0 0 0]);
    scatter3(C(1,2), C(1,3), C(1,1), [100], 'Marker', '+', 'MarkerEdgeColor', [0 0 0]);
    scatter3(featMatx(nacInds,2), featMatx(nacInds,3), featMatx(nacInds,1), 'filled', 'MarkerFaceColor', [1 0 0]);
    surf(z,y,x, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    %xlim([0 1]); ylim([0.25 1]); zlim([0 1])
    xlim([-1.5 2.5]); ylim([-1.5 2.5]); zlim([-1.5 2.5])
    %ticks = [0:0.1:1];
    %xticks(ticks); yticks(ticks); zticks(ticks);
    xlabel('z-scored spike frequency');
    ylabel('z-scored spike width');
    zlabel('z-scored latency to spike');
    suptitle('k-means clustering analysis')
    ax = gca;
    ax.CameraPosition = [25.2096  -21.9562    9.7278];
    box on;
    ax.BoxStyle = 'full';

    if labelFlag
        for i = 1:parvoInds(end)
            tmpName = [num2str(round(featMatx(i,2),2)) ', ' num2str(round(featMatx(i,3),2)) ', ' num2str(round(featMatx(i,1),2))];
            text(featMatx(i,2), featMatx(i,3), featMatx(i,1), tmpName)
        end
    end
end

normFeatTbl = table(featMatx(:,1),featMatx(:,2),featMatx(:,3), 'VariableNames', {'firstSpikeLat' 'spikeFreq' 'firstSpikeHW'}, 'RowNames',...
    [magno.traceTbl.Properties.RowNames; parvo.traceTbl.Properties.RowNames; nac.traceTbl.Properties.RowNames;]);
end

