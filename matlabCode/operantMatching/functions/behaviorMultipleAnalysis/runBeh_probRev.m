function runBeh_probRev(xlFile, animal, category)

[root, sep] = currComputer();

for i = 1:length(animal)
    [weights, dayList, ~] = xlsread(xlFile, animal{i});
    [~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
    dayList = dayList(2:end,col);
    endInd = find(cellfun(@isempty,dayList),1);
    if ~isempty(endInd)
        dayList = dayList(1:endInd-1,:);
    end

    for i = 1: length(dayList)
        sessionName = dayList{i}
        analyzeBehavioralData_operantMatching([sessionName '.asc']);
        close;
    end
end