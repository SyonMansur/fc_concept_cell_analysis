function rawTable = extractSelectiveResponses(mainFolderPath)
% parameters for zeta and windowing
max_dur = 1.5;         
resamp_num = 100;         
plot_flag = 0;                
restrict_range = [-0.5 1.0]; 
direct_quantile = false; 
jitter_size = 1;          
stitch_flag = false;         

% latency search params
radius = 0.075;  
smooth_win = 20;      

% specific cluster metrics to pull from source
cols_to_keep = {
  'csc_num', 'cluster_region', 'cluster_num', 'noiseProminence', 'firingRate', ...
  'refPeriodViolations', 'meanAmplitude', 'localMaxima', 'localMaximaRatio', ...
  'lateRangeRatio', 'meanStandardErr', 'rejectCluster', 'waveDuration', ...
  'selectivity', 'responseSelectivity'
};

rawTable = table();

% get subject folders
contents = dir(mainFolderPath);
subs = contents([contents.isdir]);
subs = subs(~ismember({subs.name}, {'.', '..'}));

for i = 1:length(subs)
    sub_folder = subs(i).name;
    % pull subject id from folder name
    sub_id_str = regexp(sub_folder, '^\d+', 'match', 'once');
    if isempty(sub_id_str), continue; end
    sub_id = str2double(sub_id_str);
    fprintf('subject: %d\n', sub_id);
    
    current_path = fullfile(mainFolderPath, sub_folder);
    % find all cluster files recursively
    mat_files = dir(fullfile(current_path, '**/clusterCharacteristics.mat'));
    
    for j = 1:length(mat_files)
        file_path = fullfile(mat_files(j).folder, mat_files(j).name);
        fprintf('  loading: %s\n', mat_files(j).name);
        
        % try to grab experiment name from path
        exp_name = '';
        try
            rel_path = erase(mat_files(j).folder, current_path);
            parts = strsplit(rel_path, filesep);
            parts = parts(~cellfun('isempty', parts));
            if ~isempty(parts), exp_name = parts{1}; end
        catch
        end
        
        clear clusterCharacteristics
        try
            load(file_path);
        catch
            continue; 
        end
        
        if ~exist('clusterCharacteristics', 'var'), continue; end
        
        file_total = 0;
        table_cols = clusterCharacteristics.Properties.VariableNames;
        
        % check each cluster in the file
        for k = 1:height(clusterCharacteristics)
            % only care if cluster shows selectivity
            if clusterCharacteristics.numSelective(k) == 0, continue; end
            
            % prep base cluster info
            base_data = table();
            for c = 1:length(cols_to_keep)
                col = cols_to_keep{c};
                if ismember(col, table_cols)
                    base_data.(col) = clusterCharacteristics.(col)(k);
                else
                    base_data.(col) = NaN; 
                end
            end
            
            % handle different field names for stimulus info
            if ismember('stim', table_cols)
                src_col = 'stim'; name_field = 'stimName';
            elseif ismember('screeningInfo', table_cols)
                src_col = 'screeningInfo'; name_field = 'name';
            else
                continue;
            end
            
            stim_array = clusterCharacteristics.(src_col){k};
            n_stim = length(stim_array);
            
            % holds info for each stimulus event
            % added spikes column so hmm can access raw data later
            stim_info = table('Size', [0 7], 'VariableTypes', {'cell', 'double', 'double', 'double', 'double', 'double', 'cell'}, ...
                             'VariableNames', {'stimName', 'originalResponseOnset', 'score', 'zetaResponseOnset', 'zetaPValue', 'num_stimuli', 'spikes'});
            
            if isempty(stim_array), continue; end
            
            for s = 1:length(stim_array)
                % check for valid spike data and metadata
                if isfield(stim_array(s), 'responseOnset') && ~isempty(stim_array(s).responseOnset) && ...
                   stim_array(s).responseOnset ~= 0 && ...
                   isfield(stim_array(s), name_field) && isfield(stim_array(s), 'score') && isfield(stim_array(s), 'spikes')
               
                    spikes = stim_array(s).spikes;
                    stitched = [];
                    cnt = -1;
                    z_p = NaN;
                    z_onset = NaN;
                    
                    % combine trials into one long spike train
                    for t = 1:length(spikes)
                        if ~isempty(spikes{t})
                            cnt = cnt + 1;
                            t_spikes = spikes{t};
                            % use one point five s window per trial
                            mask = t_spikes >= -500 & t_spikes <= 1000;
                            valid_spikes = t_spikes(mask);
                            adj_spikes = (1500 * cnt) + valid_spikes;
                            
                            if isempty(stitched)
                                stitched = adj_spikes';
                            else
                                stitched = [stitched, adj_spikes'];
                            end
                        end
                    end
                    
                    if ~isempty(stitched)
                        try
                            stitched = stitched';
                            % shift timing for seconds-based zeta input
                            onsets_ms = (0:1500:(1500 * cnt)) + 500;
                            spikes_ms = stitched + 500;
                            spikes_sec = spikes_ms / 1000;
                            onsets_sec = onsets_ms / 1000;
                            
                            % run the actual zeta test
                            [dblZetaP, sZETA, ~, sLatencies] = zetatest(spikes_sec, ...
                                onsets_sec, max_dur, resamp_num, plot_flag, ...
                                restrict_range, direct_quantile, jitter_size, stitch_flag);
                            
                            % fine tune onset using deviation derivative
                            r_onset_s = stim_array(s).responseOnset / 1000.0;
                            w_min = r_onset_s - radius;
                            w_max = r_onset_s + radius;
                            
                            % find where deviation from linear rate is steepest
                            dev = sZETA.vecRealFrac - sZETA.vecRealFracLinear;
                            s_dev = smoothdata(dev, 'movmean', smooth_win);
                            d_s = gradient(s_dev);
                            d_t = sZETA.vecSpikeT;
                            
                            v_idx = d_t >= w_min & d_t <= w_max;
                            f_onset_s = 0; 
                            
                            if any(v_idx)
                                m_slope = d_s(v_idx);
                                m_times = d_t(v_idx);
                                [~, peak_idx] = max(m_slope);
                                t_peak = m_times(peak_idx);
                                % look back for local minimum before the slope peak
                                pre_peak = find(sZETA.vecSpikeT >= w_min & sZETA.vecSpikeT <= t_peak);
                                
                                if ~isempty(pre_peak)
                                    [~, valley_idx] = min(s_dev(pre_peak));
                                    f_idx = pre_peak(valley_idx);
                                    f_onset_s = sZETA.vecSpikeT(f_idx);
                                else
                                    f_onset_s = t_peak;
                                end
                            else
                                f_onset_s = 0; 
                            end
                            
                            % use zeta default if custom search hits early floor
                            if f_onset_s <= 0.005
                                f_onset_s = sLatencies.Onset;
                            end
                            
                            z_p = dblZetaP;
                            z_onset = f_onset_s * 1000;
                       
                        catch ME
                            fprintf('zeta fail: cluster %d stim %d. %s\n', ...
                                clusterCharacteristics.cluster_num(k), s, ME.message);
                        end
                    end
                    
                    % pack stimulus result row
                    row = {stim_array(s).(name_field), stim_array(s).responseOnset, ...
                                  stim_array(s).score, z_onset, z_p, n_stim, spikes};
                    stim_info = [stim_info; row];
                end
            end
            
            if isempty(stim_info), continue; end
            file_total = file_total + height(stim_info);
            
            % label modality based on stim names
            types = cell(height(stim_info), 1);
            for r = 1:height(stim_info)
                name = lower(stim_info.stimName{r});
                if contains(name, 'audio'), types{r} = 'audio';
                elseif contains(name, 'text'), types{r} = 'text';
                else, types{r} = 'image'; end
            end
            stim_info.stimulusType = types;
            
            % merge cluster stats with specific stim results
            n_found = height(stim_info);
            base_exp = repmat(base_data, n_found, 1);
            base_exp.SubjectID = repmat(sub_id, n_found, 1);
            base_exp.Experiment = repmat({exp_name}, n_found, 1);
            
            new_rows = [base_exp, stim_info];
            
            % grow master table with outerjoin to preserve data structure
            if isempty(rawTable)
                  rawTable = new_rows;
            else
                  % standardize region column to avoid join errors
                  new_rows.cluster_region = string(new_rows.cluster_region);
                  rawTable.cluster_region = string(rawTable.cluster_region);
                  rawTable = outerjoin(rawTable, new_rows, 'MergeKeys', true);
            end
        end
        fprintf('  found %d selective responses\n', file_total);
    end
end
fprintf('extraction complete. total: %d\n', height(rawTable));
end