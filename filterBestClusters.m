function finalTable = filterBestClusters(mappedTable)
% remove zero clusters first to clean up the table
% script: filterBestClusters.m
clean = mappedTable(mappedTable.cluster_num ~= 0, :);

% grouping variables to define a unique cell/channel
groupVars = {'SubjectID', 'Experiment', 'cluster_region', 'csc_num', 'cluster_num'};

% calculate mean score for each cluster to find the most robust ones
avgScores = groupsummary(clean, groupVars, 'mean', 'score');

% identify the cluster with the highest average score per channel
% variables that define the channel itself
cscVars = {'SubjectID', 'Experiment', 'cluster_region', 'csc_num'};

% sort by score so the best cluster is at the top of each group
sorted = sortrows(avgScores, 'mean_score', 'descend');

% take the unique rows to keep only the highest scoring cluster per channel
[~, ia] = unique(sorted(:, cscVars), 'rows');
bestClusters = sorted(ia, :);

% merge back with the original mapped data to get the full response rows
finalTable = innerjoin(mappedTable, bestClusters, 'Keys', groupVars);

fprintf('filtering complete. rows: %d\n', height(finalTable));

% basic stats for the console
isFC = contains(finalTable.BroadArea, "FC");
isMTL = contains(finalTable.BroadArea, "MTL");

% print counts for unique selective units
fprintf('found %d unique selective cells in fc\n', sum(isFC));
fprintf('found %d unique selective cells in mtl\n', sum(isMTL));
end