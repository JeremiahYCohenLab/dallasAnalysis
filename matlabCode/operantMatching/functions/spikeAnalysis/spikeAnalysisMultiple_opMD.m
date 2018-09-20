function [tauList] = spikeAnalysisMultiple_opMD(xlFile)

[root, sep] = currComputer();

[~, sessionCellList, ~] = xlsread(xlFile);
sessionList = sessionCellList(2:end, 2);
cellList = sessionCellList(2:end, 1);

decayConstList = [1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192];
phasic = struct;
tonic = struct;
j = 1;
k = 1;

for i = 1:length(sessionList)
    spikeMdlTmp = spikeAnalysis_opMD(sessionList{i}, cellList{i})
    
    if any(spikeMdlTmp.(cellList{i}).postCSspikeCountRwdHist.Coefficients.pValue(2:end, 1) < 0.05)
        [phasicMinP, phasicMin] = min(spikeMdlTmp.(cellList{i}).postCSspikeCountRwdHist.Coefficients.pValue(2:end, 1));
        phasic(j).tau = decayConstList(phasicMin);
        phasic(j).session = sessionList{i};
        phasic(j).cell = cellList{i};
        j = j + 1;
    end
    if any(spikeMdlTmp.(cellList{i}).preCSspikeCountRwdHist.Coefficients.pValue(2:end, 1) < 0.05)
        [tonicMinP, tonicMin] = min(spikeMdlTmp.(cellList{i}).preCSspikeCountRwdHist.Coefficients.pValue(2:end, 1)); 
        tonic(k).tau = decayConstList(tonicMin);
        tonic(k).session = sessionList{i};
        tonic(k).cell = cellList{i};
        k = k + 1;
    end
end

tauList = struct;
tauList.phasic = phasic; 
tauList.tonic = tonic;