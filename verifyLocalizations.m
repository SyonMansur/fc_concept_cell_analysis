function finalTable = verifyLocalizations(finalTable, baseLocDir)
% establishes true location by checking clinical excel sheets
% script: verifyLocalizations.m
numRows = height(finalTable);
newColNames = {'NMM_matches', 'Anat_matches', 'aparc_aseg_matches', 'aparc_DKT_matches'};
excelColNames = {'NMM', 'Anat', 'aparc_aseg', 'aparc_DKTatlas_aseg'};

% init empty columns
emptyCellCol = cell(numRows, 1);
for col = newColNames
    finalTable.(col{1}) = emptyCellCol;
end

fprintf('processing %d rows for localization...\n', numRows);
for i = 1:numRows
    subjectID = string(finalTable.SubjectID(i));    
    clusterRegion = string(finalTable.cluster_region(i));
    
    % format target string to match electrode names in excel
    targetString = clusterRegion + "_micro";

    subFolder = "sub-" + subjectID;
    excelFile = "sub-" + subjectID + "_localizations.xlsx";
    excelFilePath = fullfile(baseLocDir, subFolder, excelFile);

    % check if clinical file exists
    if ~isfile(excelFilePath)
        for col = newColNames
            finalTable.(col{1}){i} = 'no localization';
        end
        continue;
    end
    
    try
        locTable = readtable(excelFilePath);
    catch
        for col = newColNames
            finalTable.(col{1}){i} = 'file read error';
        end
        continue;
    end
    
    % look for matching electrodes
    if ~ismember('electrode', locTable.Properties.VariableNames)
        for col = newColNames
            finalTable.(col{1}){i} = 'loc file missing electrode col';
        end
        continue;
    end
    
    electrodeCol = string(locTable.electrode);
    matchingRows = contains(electrodeCol, targetString, 'ignorecase', true);

    if ~any(matchingRows)
        for col = newColNames
            finalTable.(col{1}){i} = 'no matching electrode';
        end
    else
        % extract atlas data for matches
        for k = 1:length(newColNames)
            newCol = newColNames{k};
            excelCol = excelColNames{k};
            
            if ~ismember(excelCol, locTable.Properties.VariableNames)
                finalTable.(newCol){i} = 'loc file missing data col';
                continue;
            end

            finalTable.(newCol){i} = locTable.(excelCol)(matchingRows);
        end
    end
end

% simplify cell conflicts
% if multiple entries match, check if they agree
colsToProcess = {'NMM_matches', 'Anat_matches', 'aparc_aseg_matches', 'aparc_DKT_matches'};
for col = colsToProcess
    colName = col{1};
    for i = 1:numRows
        currentValue = finalTable.(colName){i};
        if iscell(currentValue)
            uniqueVals = unique(currentValue);
            if numel(uniqueVals) == 1
                finalTable.(colName){i} = uniqueVals{1};
            end
        end
    end
end

% map atlas keys to clean names
rawKeys = { ...
    'Left-Amygdala', 'Left-Cerebral-White-Matter', 'Left-Hippocampus', ...
    'Left-Lateral-Ventricle', 'Left-Putamen', ...
    'Right-Amygdala', 'Right-Cerebral-White-Matter', 'Right-Hippocampus', ...
    'Right-Pallidum', 'Right-Putamen', ...
    'Unknown', 'unknown', ... 
    'ctx-lh-caudalanteriorcingulate', 'ctx-lh-cuneus', 'ctx-lh-entorhinal', ...
    'ctx-lh-fusiform', 'ctx-lh-inferiortemporal', 'ctx-lh-insula', ...
    'ctx-lh-isthmuscingulate', 'ctx-lh-lateralorbitofrontal', 'ctx-lh-lingual', ...
    'ctx-lh-medialorbitofrontal', 'ctx-lh-parahippocampal', 'ctx-lh-parsopercularis', ...
    'ctx-lh-parstriangularis', 'ctx-lh-precuneus', 'ctx-lh-rostralanteriorcingulate', ...
    'ctx-lh-superiorfrontal', 'ctx-lh-superiortemporal', 'ctx-lh-transversetemporal', ...
    'ctx-rh-caudalanteriorcingulate', 'ctx-rh-entorhinal', 'ctx-rh-fusiform', ...
    'ctx-rh-insula', 'ctx-rh-isthmuscingulate', 'ctx-rh-lingual', ...
    'ctx-rh-medialorbitofrontal', 'ctx-rh-parahippocampal', 'ctx-rh-parsopercularis', ...
    'ctx-rh-parstriangularis', 'ctx-rh-pericalcarine', 'ctx-rh-posteriorcingulate', ...
    'ctx-rh-precuneus', 'ctx-rh-rostralanteriorcingulate', 'ctx-rh-superiorfrontal', ...
    'ctx-rh-superiortemporal', 'ctx-rh-transversetemporal' ...
};

cleanValues = rawKeys;
for k = 1:length(rawKeys)
    val = erase(rawKeys{k}, {'Left-', 'Right-', 'ctx-lh-', 'ctx-rh-'});
    val = replace(val, '-', ' ');
    if ~isempty(val)
        val(1) = upper(val(1));
    end
    cleanValues{k} = val;
end
locMap = containers.Map(rawKeys, cleanValues);

% resolve trueLocation
finalTable.trueLocation = strings(numRows, 1);
finalTable.verified = false(numRows, 1);

for i = 1:numRows
    locCandidate = finalTable.aparc_DKT_matches{i};
    useLocalization = false;
    finalLocName = "";
    
    if ischar(locCandidate) || isstring(locCandidate)
        strCand = string(locCandidate);
        isErrorMsg = contains(strCand, ["no localization", "file read error", "missing"], 'ignorecase', true);
        
        if ~isErrorMsg && isKey(locMap, strCand)
            mappedName = locMap(strCand);
            if ~strcmpi(mappedName, 'Unknown')
                useLocalization = true;
                finalLocName = mappedName;
            end
        end
    end
    
    if useLocalization
        finalTable.trueLocation(i) = finalLocName;
        finalTable.verified(i) = true;
    else
        % fallback to mapped region if clinical loc is missing
        if ismember('CleanedRegion', finalTable.Properties.VariableNames) && ~isempty(finalTable.CleanedRegion(i))
            finalTable.trueLocation(i) = string(finalTable.CleanedRegion(i));
        else
            finalTable.trueLocation(i) = "Unknown";
        end
    end
end

fprintf('localization complete. verified: %d | unverified: %d\n', ...
    sum(finalTable.verified), sum(~finalTable.verified));
end