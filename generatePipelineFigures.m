function generatePipelineFigures(finalTable)
% =========================================================================
% all of the figures!!
% adapted for pipeline usage 12/13/25
% =========================================================================

if nargin < 1 || isempty(finalTable)
    error('finalTable is missing. check your pipeline flow.');
end

zeta_data = finalTable;

score_threshold = 4.25;
time_limits = [0 1000];

% round latencies
if isnumeric(zeta_data.zetaResponseOnset)
    zeta_data.zetaResponseOnset = round(zeta_data.zetaResponseOnset);
end

raw_mask = zeta_data.zetaResponseOnset > time_limits(1) & ...
           zeta_data.zetaResponseOnset <= time_limits(2);
clean_data = zeta_data(raw_mask, :);

raw_FC  = clean_data(clean_data.BroadArea == "FC", :);
raw_MTL = clean_data(clean_data.BroadArea == "MTL", :);

% using truelocation from verification step
raw_HC  = clean_data(contains(clean_data.trueLocation, "Hippocampus", "IgnoreCase", true), :);
raw_AMY = clean_data(contains(clean_data.trueLocation, "amygdala", "IgnoreCase", true), :);
raw_EC  = clean_data(contains(clean_data.trueLocation, "entorhi", "IgnoreCase", true), :);
raw_PHC = clean_data(contains(clean_data.trueLocation, "parahippo", "IgnoreCase", true), :);

raw_ACC  = clean_data(contains(clean_data.trueLocation, "anterior cingulate", "IgnoreCase", true), :);
raw_OFC  = clean_data(contains(clean_data.trueLocation, "orbitofrontal", "IgnoreCase", true), :);
raw_OPER = clean_data(contains(clean_data.trueLocation, "operculum", "IgnoreCase", true), :);
raw_MOT  = clean_data(contains(clean_data.trueLocation, "motor", "IgnoreCase", true), :);

cell_FC  = aggregate_cells(raw_FC);
cell_MTL = aggregate_cells(raw_MTL);
cell_HC  = aggregate_cells(raw_HC);
cell_AMY = aggregate_cells(raw_AMY);
cell_EC  = aggregate_cells(raw_EC);
cell_PHC = aggregate_cells(raw_PHC);
cell_ACC = aggregate_cells(raw_ACC);
cell_OFC = aggregate_cells(raw_OFC);
cell_OPER= aggregate_cells(raw_OPER);
cell_MOT = aggregate_cells(raw_MOT);

% =========================================================================
% groups

mtl_groups = {
   cell_MTL, raw_MTL, 'MTL (General)',   [0.00 0.57 0.55];
   cell_PHC, raw_PHC, 'PHC',             [0.49 0.18 0.56]; 
   cell_EC,  raw_EC,  'EC',              [0.47 0.67 0.19]; 
   cell_HC,  raw_HC,  'Hipp',            [0.00 0.45 0.74]; 
   cell_AMY, raw_AMY, 'Amy',             [0.85 0.33 0.10]; 
};

fc_groups = {
   cell_FC,   raw_FC,   'FC (General)',       [0.63 0.07 0.18];
   cell_MOT,  raw_MOT,  'Motor Area',         [0.30 0.74 0.93]; 
   cell_ACC,  raw_ACC,  'Anterior Cingulate', [0.92 0.69 0.12];
   cell_OFC,  raw_OFC,  'Orbitofrontal',      [0.49 0.18 0.56];
   cell_OPER, raw_OPER, 'Operculum',          [0.47 0.67 0.19];
};

% =========================================================================
% =========================================================================
fprintf('Generating Figures...\n');

plot_histograms(mtl_groups, 'MTL Latency (Cell Median)', 'cell');
plot_histograms(fc_groups,  'FC Latency (Cell Median)',  'cell');
plot_histograms(mtl_groups, 'MTL Latency (All Responses)', 'raw');
plot_histograms(fc_groups,  'FC Latency (All Responses)', 'raw');

mtl_sub = mtl_groups(2:end, :);
fc_sub  = fc_groups(2:end, :); 
plot_bars(mtl_sub, 'Mean Latency by MTL Region', 'latency');
plot_bars(fc_sub,  'Mean Latency by FC Region',  'latency');
plot_bars(mtl_sub, 'Mean Responsiveness by MTL Region', 'count');
plot_bars(fc_sub,  'Mean Responsiveness by FC Region',  'count');

plot_scatter_zoomed(mtl_groups, 'MTL Selectivity vs Latency');
plot_scatter_zoomed(fc_groups,  'FC Selectivity vs Latency');

plot_cdfs(mtl_sub, 'MTL CDF Comparison');
plot_cdfs(fc_sub,  'FC CDF Comparison');

% figures combined
plot_combined_hist_bar(mtl_groups, 'MTL Latency by Region', 'latency');
plot_combined_hist_bar(fc_groups,  'FC Latency by Region',  'latency');
plot_combined_hist_bar(mtl_groups, 'MTL Responsiveness by Region', 'count');
plot_combined_hist_bar(fc_groups,  'FC Responsiveness by Region',  'count');
plot_combined_hist_bar(mtl_groups, 'MTL Proportion by Region', 'proportion');
plot_combined_hist_bar(fc_groups,  'FC Proportion by Region',  'proportion');

% =========================================================================
% fc split
% =========================================================================
fprintf('Running Hypothesis Test...\n');
if ~isempty(raw_FC) && ~isempty(raw_MTL)
   strong_fc = raw_FC.zetaResponseOnset(raw_FC.score >= score_threshold);
   weak_fc   = raw_FC.zetaResponseOnset(raw_FC.score < score_threshold);
   mtl_all   = raw_MTL.zetaResponseOnset;
   n_strong = length(strong_fc); n_weak = length(weak_fc); n_mtl = length(mtl_all);
   
   figure('Name', 'Hypothesis: Tri-Split', 'Color', 'w', 'Position', [100, 100, 900, 500]);
   hold on;
   histogram(strong_fc, 'BinWidth', 40, 'FaceColor', [0.8 0 0], 'FaceAlpha', 0.6, 'EdgeColor', 'k');
   histogram(mtl_all,   'BinWidth', 40, 'FaceColor', [0 0.6 0.6], 'FaceAlpha', 0.4, 'EdgeColor', 'k');
   histogram(weak_fc,   'BinWidth', 40, 'FaceColor', [0 0 0.6], 'FaceAlpha', 0.4, 'EdgeColor', 'k');
   if ~isempty(strong_fc), xline(median(strong_fc), '--r', 'LineWidth', 2); end
   if ~isempty(mtl_all),   xline(median(mtl_all), '--c', 'LineWidth', 2); end
   if ~isempty(weak_fc),   xline(median(weak_fc),  '--b', 'LineWidth', 2); end
   legend({sprintf('Strong FC (N=%d)', n_strong), sprintf('MTL (N=%d)', n_mtl), sprintf('Weak FC (N=%d)', n_weak)});
   title(sprintf('FC Population Split vs MTL (Score Thresh = %.1f)', score_threshold));
   xlabel('Response Latency (ms)'); ylabel('Number of Responses');
   hold off;
   
   figure('Name', 'Hypothesis: Strong vs Weak FC', 'Color', 'w', 'Position', [150, 150, 900, 500]);
   hold on;
   histogram(strong_fc, 'BinWidth', 40, 'FaceColor', [0.8 0 0], 'FaceAlpha', 0.6, 'EdgeColor', 'k');
   histogram(weak_fc,   'BinWidth', 40, 'FaceColor', [0 0 0.6], 'FaceAlpha', 0.4, 'EdgeColor', 'k');
   if ~isempty(strong_fc), xline(median(strong_fc), '--r', 'LineWidth', 2); end
   if ~isempty(weak_fc),   xline(median(weak_fc),  '--b', 'LineWidth', 2); end
   [~, p_val_sw] = ttest2(strong_fc, weak_fc);
   legend({sprintf('Strong FC (N=%d)', n_strong), sprintf('Weak FC (N=%d)', n_weak)});
   title(sprintf('Strong FC vs Weak FC (p = %.1e)', p_val_sw));
   xlabel('Response Latency (ms)'); ylabel('Number of Responses');
   hold off;
   
   figure('Name', 'Hypothesis: Sorted Boxplot', 'Color', 'w', 'Position', [100, 100, 700, 500]);
   m_strong = mean(strong_fc); m_weak = mean(weak_fc); m_mtl = mean(mtl_all);
   data_cell = {strong_fc, mtl_all, weak_fc};
   names_cell = {'Strong FC', 'MTL', 'Weak FC'};
   means_arr  = [m_strong, m_mtl, m_weak];
   [sorted_means, sort_idx] = sort(means_arr, 'ascend');
   sorted_data = data_cell(sort_idx); sorted_names = names_cell(sort_idx);
   all_plot_data = []; all_labels = [];
   for i = 1:3
       curr_data = sorted_data{i};
       all_plot_data = [all_plot_data; curr_data];
       all_labels = [all_labels; repmat(categorical(sorted_names(i)), length(curr_data), 1)];
   end
   all_labels = reordercats(all_labels, sorted_names);
   boxchart(all_labels, all_plot_data, 'BoxFaceColor', [0.3 0.3 0.3], 'MarkerStyle', '.'); hold on;
   plot(1:3, sorted_means, 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 12, 'LineWidth', 2);
   ylabel('Response Time (ms)'); title('Latency by Group (Ordered by Mean)');
   if ~isempty(all_plot_data)
       y_limit_percentile = prctile(all_plot_data, 98); y_max = y_limit_percentile * 1.2; ylim([0 y_max]);
       pval_step = y_max * 0.08; y_start_pval = y_max * 0.85;
       [~, p1] = ttest2(sorted_data{1}, sorted_data{2});
       line([1 2], [y_start_pval y_start_pval], 'Color', 'k');
       text(1.5, y_start_pval + pval_step*0.3, sprintf('p=%.1e', p1), 'HorizontalAlignment', 'center');
       [~, p2] = ttest2(sorted_data{2}, sorted_data{3});
       y_line2 = y_start_pval + pval_step;
       line([2 3], [y_line2 y_line2], 'Color', 'k');
       text(2.5, y_line2 + pval_step*0.3, sprintf('p=%.1e', p2), 'HorizontalAlignment', 'center');
   end
   hold off;
end

% =========================================================================
% scatters
% =========================================================================
plot_general_stats(cell_MTL, cell_FC);
plot_score_scatter(cell_MTL, cell_FC);
plot_raw_score_scatter(raw_MTL, raw_FC);

% =========================================================================
% paired
% =========================================================================
fprintf('\n--- Running Paired Analysis ---\n');
[G, subjects, stims] = findgroups(clean_data.SubjectID, clean_data.stimName);
has_mtl = splitapply(@(x) any(x == "MTL"), clean_data.BroadArea, G);
has_fc  = splitapply(@(x) any(x == "FC"), clean_data.BroadArea, G);
valid_pairs = has_mtl & has_fc;
valid_indices = find(valid_pairs);
num_pairs = length(valid_indices);
fprintf('Found %d paired Subject-Stimulus interactions.\n', num_pairs);

if num_pairs > 0
   med_mtl_onsets = nan(num_pairs, 1); med_fc_onsets = nan(num_pairs, 1);
   med_fc_scores = nan(num_pairs, 1); med_mtl_scores = nan(num_pairs, 1);
   
   for k = 1:num_pairs
       idx = valid_indices(k);
       curr_mask = (G == idx);
       curr_table = clean_data(curr_mask, :);
       mtl_sub = curr_table(curr_table.BroadArea == "MTL", :);
       fc_sub  = curr_table(curr_table.BroadArea == "FC", :);
       med_mtl_onsets(k) = median(mtl_sub.zetaResponseOnset, 'omitnan');
       med_fc_onsets(k)  = median(fc_sub.zetaResponseOnset, 'omitnan');
       med_fc_scores(k)  = median(fc_sub.score, 'omitnan');
       med_mtl_scores(k) = median(mtl_sub.score, 'omitnan');
   end
   
   paired_tbl = table(med_mtl_onsets, med_fc_onsets, med_fc_scores, med_mtl_scores);
   paired_tbl.Diff = paired_tbl.med_fc_onsets - paired_tbl.med_mtl_onsets;
   
   figure('Name', 'Paired Analysis: Latency', 'Color', 'w', 'Position', [100, 100, 600, 500]);
   box_data = [paired_tbl.med_mtl_onsets, paired_tbl.med_fc_onsets];
   boxplot(box_data, 'Labels', {'MTL', 'FC'}, 'Symbol', ''); hold on;
   x_swarm = [ones(num_pairs, 1); 2*ones(num_pairs, 1)];
   y_swarm = [paired_tbl.med_mtl_onsets; paired_tbl.med_fc_onsets];
   swarmchart(x_swarm, y_swarm, 15, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.5);
   ylabel('Median Response Onset (ms)');
   [~, p_paired] = ttest(paired_tbl.med_mtl_onsets, paired_tbl.med_fc_onsets);
   title({sprintf('Paired Median Onsets (N=%d Pairs)', num_pairs), sprintf('p=%.1e', p_paired)});
   hold off;
   
   figure('Name', 'Paired Analysis: Difference', 'Color', 'w', 'Position', [700, 100, 400, 500]);
   boxplot(paired_tbl.Diff, 'Labels', {'FC - MTL'}, 'Symbol', '.'); yline(0, '--k');
   [~, p_diff] = ttest(paired_tbl.Diff);
   title({'Distribution of Pairwise Differences', sprintf('Mean Diff vs 0: p = %.1e', p_diff)});
   
   % split by Score
   is_high_fc = paired_tbl.med_fc_scores >= score_threshold;
   is_low_fc  = paired_tbl.med_fc_scores < score_threshold;
   high_fc_pairs = paired_tbl(is_high_fc, :); low_fc_pairs = paired_tbl(is_low_fc, :);
   
   figure('Name', 'Paired Split by FC Score', 'Color', 'w', 'Position', [100, 100, 1000, 500]);
   t = tiledlayout(1, 2, 'Padding', 'compact');
   nexttile;
   histogram(high_fc_pairs.med_fc_onsets, 'BinWidth', 50, 'FaceColor', [0.8 0 0], 'FaceAlpha', 0.6); hold on;
   histogram(low_fc_pairs.med_fc_onsets, 'BinWidth', 50, 'FaceColor', [0 0 0.6], 'FaceAlpha', 0.6);
   xline(median(high_fc_pairs.med_fc_onsets), '--r'); xline(median(low_fc_pairs.med_fc_onsets), '--b');
   legend({'Strong FC Pairs', 'Weak FC Pairs'}); title('FC Onset: Split by Score'); hold off;
   
   nexttile;
   histogram(high_fc_pairs.med_mtl_onsets, 'BinWidth', 50, 'FaceColor', [0.8 0 0], 'FaceAlpha', 0.6); hold on;
   histogram(low_fc_pairs.med_mtl_onsets, 'BinWidth', 50, 'FaceColor', [0 0 0.6], 'FaceAlpha', 0.6);
   xline(median(high_fc_pairs.med_mtl_onsets), '--r'); xline(median(low_fc_pairs.med_mtl_onsets), '--b');
   legend({'MTL (w/ Strong FC)', 'MTL (w/ Weak FC)'}); title('MTL Onset: Split by FC Partner Score'); hold off;

else
   fprintf('No valid pairs found.\n');
end

max_thresh = 10;
plot_median_vs_threshold(raw_MTL, 'MTL General', max_thresh, [0.00 0.57 0.55]);
plot_median_vs_threshold(raw_FC,  'FC General',  max_thresh, [0.63 0.07 0.18]);
plot_broad_area_pie(clean_data);

% =========================================================================
% paired latency
% =========================================================================
if exist('paired_tbl', 'var') && ~isempty(paired_tbl)
    figure('Name', 'Paired Latency Correlation', 'Color', 'w', 'Position', [100, 100, 600, 600]);
    scatter(paired_tbl.med_mtl_onsets, paired_tbl.med_fc_onsets, 25, 'filled', ...
            'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerFaceAlpha', 0.4); hold on;
    max_val = max([paired_tbl.med_mtl_onsets; paired_tbl.med_fc_onsets]);
    plot([0 max_val], [0 max_val], '--k', 'LineWidth', 1.5);
    l = lsline; l.Color = 'r'; l.LineWidth = 2;
    [rho, p_val] = corr(paired_tbl.med_mtl_onsets, paired_tbl.med_fc_onsets, 'Type', 'Spearman');
    title({'MTL vs FC Onset Synchronization', sprintf('Spearman r = %.2f, p = %.2e', rho, p_val)});
    xlabel('MTL Onset (ms)'); ylabel('FC Onset (ms)'); axis square;
end

% =========================================================================
% amplitude coupling
% =========================================================================
if exist('paired_tbl', 'var') && ~isempty(paired_tbl)
    figure('Name', 'Amplitude Coupling', 'Color', 'w', 'Position', [700, 100, 600, 600]);
    
    scatter(paired_tbl.med_mtl_scores, paired_tbl.med_fc_scores, 25, 'filled', ...
        'MarkerFaceColor', [0.5 0 0.5], 'MarkerFaceAlpha', 0.4); 
    hold on;
    xline(score_threshold, '--k'); yline(score_threshold, '--k');
    
    
    l = lsline; l.Color = 'r'; l.LineWidth = 2;
    
    xlabel('MTL Score'); ylabel('FC Score');
    title('Response Strength Coupling');
    xlim([0 10]); ylim([0 10]);
    
    [r_amp, p_amp] = corr(paired_tbl.med_mtl_scores, paired_tbl.med_fc_scores, 'Type', 'Spearman');
    subtitle(sprintf('Correlation: r = %.2f, p = %.2e', r_amp, p_amp));
    % text(8, 9, 'Co-Active', 'Color', 'r', 'FontWeight', 'bold');
    % text(8, 1, 'MTL Only', 'Color', 'b', 'FontWeight', 'bold');
    grid on; box off; axis square;
end

% =========================================================================
% lead lag
% =========================================================================
if exist('paired_tbl', 'var') && ~isempty(paired_tbl)
    figure('Name', 'Lead-Lag Distribution', 'Color', 'w', 'Position', [100, 100, 800, 400]);
    deltas = paired_tbl.med_fc_onsets - paired_tbl.med_mtl_onsets;
    histogram(deltas, 'BinWidth', 20, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'k'); hold on;
    xline(0, 'r', 'LineWidth', 2);
    mtl_leads_pct = sum(deltas > 0) / length(deltas) * 100;
    fc_leads_pct = sum(deltas < 0) / length(deltas) * 100;
    title('Distribution of Relative Timing (FC Onset - MTL Onset)');
    stats_txt = {sprintf('MTL Leads: %.1f%%', mtl_leads_pct), sprintf('Median Lag: %.1f ms', median(deltas))};
    annotation('textbox', [0.15 0.7 0.2 0.2], 'String', stats_txt, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
end

% =========================================================================
% subregion pairwise
% =========================================================================
fprintf('Generating Sub-Region Matrix...\n');
mtl_subs = {'Hippocampus', 'Amygdala', 'Entorhinal', 'Parahippocampal'};
fc_subs  = {'Anterior Cingulate', 'Orbitofrontal', 'Motor', 'Operculum'};
res_matrix = nan(length(mtl_subs), length(fc_subs));
for m = 1:length(mtl_subs)
    for f = 1:length(fc_subs)
        idx_m = contains(clean_data.CleanedRegion, mtl_subs{m}, 'IgnoreCase', true);
        idx_f = contains(clean_data.CleanedRegion, fc_subs{f},  'IgnoreCase', true);
        sub_data = clean_data(idx_m | idx_f, :);
        [G_sub, ~, ~] = findgroups(sub_data.SubjectID, sub_data.stimName);
        has_m_sub = splitapply(@(x) any(contains(x, mtl_subs{m}, 'IgnoreCase', true)), sub_data.CleanedRegion, G_sub);
        has_f_sub = splitapply(@(x) any(contains(x, fc_subs{f},  'IgnoreCase', true)), sub_data.CleanedRegion, G_sub);
        valid_sub_pairs = has_m_sub & has_f_sub;
        if sum(valid_sub_pairs) > 5 
            pair_indices = find(valid_sub_pairs); diffs = [];
            for k = 1:length(pair_indices)
                curr_mask = (G_sub == pair_indices(k)); curr_t = sub_data(curr_mask, :);
                lat_m = median(curr_t.zetaResponseOnset(contains(curr_t.CleanedRegion, mtl_subs{m},'IgnoreCase',true)));
                lat_f = median(curr_t.zetaResponseOnset(contains(curr_t.CleanedRegion, fc_subs{f},'IgnoreCase',true)));
                diffs(end+1) = lat_f - lat_m;
            end
            res_matrix(m, f) = mean(diffs); 
        end
    end
end
figure('Name', 'Sub-Region Connectivity Lag', 'Color', 'w');
h = heatmap(fc_subs, mtl_subs, res_matrix); h.Colormap = jet; h.MissingDataColor = [0.9 0.9 0.9];
h.Title = 'Mean Latency Lag (FC - MTL) [ms]';

% =========================================================================
% best response
% =========================================================================
fprintf('Running Detailed Best Response Analysis...\n');

% index
clean_data.OrigIdx = (1:height(clean_data))';

% group by subject / region
[G_best, ~] = findgroups(clean_data.SubjectID, clean_data.CleanedRegion);
best_idxs = splitapply(@(s, idx) idx(find(s == max(s), 1)), clean_data.score, clean_data.OrigIdx, G_best);
best_table = clean_data(best_idxs, :);

% medianlat column
best_table.MedianLat = best_table.zetaResponseOnset;

% MTL
best_MTL = best_table(best_table.BroadArea == "MTL", :);
best_HC  = best_table(contains(best_table.CleanedRegion, "Hippocampus", "IgnoreCase", true), :);
best_AMY = best_table(contains(best_table.CleanedRegion, "amygdala", "IgnoreCase", true), :);
best_EC  = best_table(contains(best_table.CleanedRegion, "entorhi", "IgnoreCase", true), :);
best_PHC = best_table(contains(best_table.CleanedRegion, "parahippo", "IgnoreCase", true), :);
% FC
best_FC  = best_table(best_table.BroadArea == "FC", :);
best_ACC = best_table(contains(best_table.CleanedRegion, "anterior cingulate", "IgnoreCase", true), :);
best_OFC = best_table(contains(best_table.CleanedRegion, "orbitofrontal", "IgnoreCase", true), :);
best_OPER= best_table(contains(best_table.CleanedRegion, "operculum", "IgnoreCase", true), :);
best_MOT = best_table(contains(best_table.CleanedRegion, "motor", "IgnoreCase", true), :);

% plotting groups
best_mtl_groups = {
   best_MTL, [], 'MTL (Best)',   [0.00 0.57 0.55];
   best_PHC, [], 'PHC (Best)',   [0.49 0.18 0.56];
   best_EC,  [], 'EC (Best)',    [0.47 0.67 0.19];
   best_HC,  [], 'Hipp (Best)',  [0.00 0.45 0.74];
   best_AMY, [], 'Amy (Best)',   [0.85 0.33 0.10];
};

best_fc_groups = {
   best_FC,   [], 'FC (Best)',       [0.63 0.07 0.18];
   best_MOT,  [], 'Motor (Best)',    [0.30 0.74 0.93];
   best_ACC,  [], 'ACC (Best)',      [0.92 0.69 0.12];
   best_OFC,  [], 'OFC (Best)',      [0.49 0.18 0.56];
   best_OPER, [], 'Oper (Best)',     [0.47 0.67 0.19];
};

plot_combined_hist_bar(best_mtl_groups, 'Best Response Latencies: MTL Regions', 'latency');
plot_combined_hist_bar(best_fc_groups,  'Best Response Latencies: FC Regions',  'latency');

end

% =========================================================================
% helper functions
% =========================================================================

function [uTable] = aggregate_cells(t)
   if isempty(t), uTable = table(); return; end
   vars = {'SubjectID', 'Experiment', 'cluster_region', 'csc_num', 'CleanedRegion', 'num_stimuli'};
   valid_rows = true(height(t), 1);
   for v = 1:length(vars), col = t.(vars{v}); if isnumeric(col), valid_rows = valid_rows & ~isnan(col); end; end
   t = t(valid_rows, :); if isempty(t), uTable = table(); return; end
   [g, uTable] = findgroups(t(:, vars));
   uTable.NumResponses = splitapply(@numel, t.score, g);
   uTable.MedianLat    = splitapply(@(x) median(x,'omitnan'), t.zetaResponseOnset, g);
   uTable.MeanScore    = splitapply(@mean, t.score, g);
   uTable.BaseFR       = splitapply(@mean, t.firingRate, g);
   uTable.num_stimuli  = splitapply(@mean, t.num_stimuli, g);
end

function plot_broad_area_pie(tbl)
   figure('Name', 'Broad Area Distribution', 'Color', 'w', 'Position', [100, 100, 600, 500]);
   vars = {'SubjectID', 'Experiment', 'cluster_region', 'csc_num', 'CleanedRegion'};
   [~, uniqueIdx] = unique(tbl(:, vars), 'rows'); uniqueCells = tbl(uniqueIdx, :);
   areas = uniqueCells.BroadArea; [groups, ~, idx] = unique(areas); counts = accumarray(idx, 1);
   [counts, sortIdx] = sort(counts, 'descend'); groups = groups(sortIdx);
   colors = zeros(length(groups), 3); base_cmap = lines(length(groups));
   for i = 1:length(groups)
       if strcmpi(groups(i), 'MTL'), colors(i,:) = [0.00 0.57 0.55]; 
       elseif strcmpi(groups(i), 'FC'), colors(i,:) = [0.63 0.07 0.18]; 
       else, colors(i,:) = base_cmap(i, :); end
   end
   p = pie(counts);
   for k = 2:2:length(p)
       label_idx = k/2; pct = counts(label_idx) / sum(counts) * 100;
       if pct > 1, p(k).String = sprintf('%s\n(N=%d, %.1f%%)', groups(label_idx), counts(label_idx), pct); p(k).FontSize = 10; p(k).FontWeight = 'bold';
       else, p(k).String = ''; end
   end
   for k = 1:2:length(p), patch_idx = (k+1)/2; p(k).FaceColor = colors(patch_idx, :); p(k).EdgeColor = 'k'; end
   title(sprintf('Distribution of Recorded Cells (Total N=%d)', sum(counts)), 'FontSize', 14);
end

function plot_general_stats(mtl, fc)
   figure('Name', 'MTL vs FC General Stats', 'Color', 'w', 'Position', [50, 50, 1400, 400]);
   t = tiledlayout(1, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
   metrics = {'MedianLat', 'NumResponses', 'Proportion', 'BaseFR'}; titles = {'Latency (ms)', 'Responsiveness (Count)', 'Responsiveness (Prop)', 'Baseline FR (Hz)'};
   for i = 1:4
       nexttile;
       if strcmp(metrics{i}, 'Proportion'), d_mtl = mtl.NumResponses ./ mtl.num_stimuli; d_fc = fc.NumResponses ./ fc.num_stimuli;
       else, d_mtl = mtl.(metrics{i}); d_fc = fc.(metrics{i}); end
       d_mtl = d_mtl(~isnan(d_mtl)); d_fc = d_fc(~isnan(d_fc));
       all_data = [d_mtl; d_fc]; labels = [repmat(categorical({'MTL'}), length(d_mtl), 1); repmat(categorical({'FC'}), length(d_fc), 1)];
       boxchart(labels, all_data, 'BoxFaceColor', [0.3 0.3 0.3], 'MarkerStyle', '.'); hold on;
       m_mtl = mean(d_mtl); m_fc = mean(d_fc); plot(1, m_mtl, 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.5);
       plot(2, m_fc, 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.5);
       [~, p] = ttest2(d_mtl, d_fc);
       if ~isempty(all_data), y_max = prctile(all_data, 98) * 1.2; ylim([0 y_max]); line([1 2], [y_max*0.85 y_max*0.85], 'Color', 'k'); text(1.5, y_max*0.92, sprintf('p=%.1e', p), 'HorizontalAlignment', 'center', 'FontSize', 10); end
       title(titles{i}); box off; hold off;
   end
   title(t, 'MTL vs FC General Comparison', 'FontSize', 14, 'FontWeight', 'bold');
end

function plot_score_scatter(mtl, fc)
   figure('Name', 'Score vs Latency (General)', 'Color', 'w', 'Position', [50, 550, 800, 400]);
   t = tiledlayout(1, 2, 'Padding', 'compact');
   nexttile; scatter(mtl.MedianLat, mtl.MeanScore, 20, 'filled', 'MarkerFaceColor', [0.00 0.57 0.55], 'MarkerFaceAlpha', 0.6);
   hold on; lsline; [r, p] = corr(mtl.MedianLat, mtl.MeanScore, 'Type', 'Spearman');
   title(sprintf('MTL General (r=%.2f, p=%.1e)', r, p)); xlabel('Latency (ms)'); ylabel('Mean Score'); xlim([0 1000]); box off; hold off;
   nexttile; scatter(fc.MedianLat, fc.MeanScore, 20, 'filled', 'MarkerFaceColor', [0.63 0.07 0.18], 'MarkerFaceAlpha', 0.6);
   hold on; lsline; [r, p] = corr(fc.MedianLat, fc.MeanScore, 'Type', 'Spearman');
   title(sprintf('FC General (r=%.2f, p=%.1e)', r, p)); xlabel('Latency (ms)'); ylabel('Mean Score'); xlim([0 1000]); box off; hold off;
end

function plot_raw_score_scatter(raw_mtl, raw_fc)
   figure('Name', 'All Responses: Score vs Latency', 'Color', 'w', 'Position', [100, 100, 800, 400]);
   t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
   nexttile; lat = raw_mtl.zetaResponseOnset; sc = raw_mtl.score;
   scatter(lat, sc, 10, 'filled', 'MarkerFaceColor', [0.00 0.57 0.55], 'MarkerFaceAlpha', 0.3); hold on; l = lsline; l.Color = 'k'; l.LineWidth = 1.5;
   [r, p] = corr(lat, sc, 'Type', 'Spearman'); title(sprintf('MTL (All Responses)\nr=%.2f, p=%.1e', r, p)); xlabel('Latency (ms)'); ylabel('ZETA Score'); xlim([0 1000]); box off; grid off; set(gca, 'TickDir', 'out'); hold off;
   nexttile; lat = raw_fc.zetaResponseOnset; sc = raw_fc.score;
   scatter(lat, sc, 10, 'filled', 'MarkerFaceColor', [0.63 0.07 0.18], 'MarkerFaceAlpha', 0.3); hold on; l = lsline; l.Color = 'k'; l.LineWidth = 1.5;
   [r, p] = corr(lat, sc, 'Type', 'Spearman'); title(sprintf('FC (All Responses)\nr=%.2f, p=%.1e', r, p)); xlabel('Latency (ms)'); ylabel('ZETA Score'); xlim([0 1000]); box off; grid off; set(gca, 'TickDir', 'out'); hold off;
   title(t, 'Score vs Latency (All Response Events)', 'FontSize', 14, 'FontWeight', 'bold');
end

function plot_median_vs_threshold(tbl, regionName, maxThresh, col)
   thresholds = 3:1:maxThresh; medians = nan(size(thresholds)); counts = nan(size(thresholds));
   for i = 1:length(thresholds)
       thr = thresholds(i); subset = tbl.zetaResponseOnset(tbl.score >= thr);
       if ~isempty(subset), medians(i) = median(subset); counts(i) = length(subset); end
   end
   figure('Name', ['Median vs Threshold: ' regionName], 'Color', 'w');
   plot(thresholds, medians, '-o', 'LineWidth', 2, 'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
   title(sprintf('%s: Median Latency vs Score Threshold', regionName), 'FontSize', 14); xlabel('Score Threshold (>=)', 'FontSize', 12); ylabel('Median Response Time (ms)', 'FontSize', 12); grid on; box off;
   for i = 1:length(thresholds), if ~isnan(medians(i)), text(thresholds(i), medians(i) + 15, sprintf('n=%d', counts(i)), 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'k'); end; end
   xlim([2.5 maxThresh+0.5]); ylim([min(medians)*0.8 max(medians)*1.2]); set(gca, 'TickDir', 'out', 'FontSize', 12);
end

function plot_histograms(groups, figTitle, mode)
   figure('Name', figTitle, 'Color', 'w', 'Position', [50, 50, 1200, 400]);
   count = size(groups, 1); t = tiledlayout(1, count, 'Padding', 'compact', 'TileSpacing', 'compact');
   for i = 1:count
       nexttile; col = groups{i, 4}; name = groups{i, 3};
       if strcmp(mode, 'cell'), tbl = groups{i, 1}; data = tbl.MedianLat; lbl='Median Latency (ms)'; ylbl='Number of Cells';
       else, tbl = groups{i, 2}; data = tbl.zetaResponseOnset; lbl='Response Time (ms)'; ylbl='Number of Responses'; end
       if isempty(data), title([name ' (No Data)']); continue; end
       histogram(data, 'BinWidth', 50, 'FaceColor', col, 'FaceAlpha', 0.7, 'EdgeColor', 'k');
       med_val = median(data); mean_val = mean(data); xline(med_val, '--k', 'LineWidth', 1.5);
       stats_str = sprintf('N=%d\nMed: %.0f\nMean: %.0f', length(data), med_val, mean_val);
       text(0.95, 0.95, stats_str, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'none', 'Margin', 1);
       title(name, 'FontSize', 11); xlabel(lbl); if i==1, ylabel(ylbl); end
       xlim([0 1000]); box on;
   end
   title(t, figTitle, 'FontSize', 16);
end

function plot_bars(groups, figTitle, mode)
   figure('Name', figTitle, 'Color', 'w', 'Position', [100, 100, 800, 500]); hold on;
   n_groups = size(groups, 1); means = zeros(1, n_groups); sems = zeros(1, n_groups);
   for i = 1:n_groups
       tbl = groups{i, 1}; if isempty(tbl), continue; end
       if strcmp(mode, 'latency'), d = tbl.MedianLat; lbl='Response Time (ms)'; else, d = tbl.NumResponses; lbl='Avg Responses / Cell'; end
       if isempty(d), continue; end
       means(i) = mean(d); sems(i) = std(d)/sqrt(length(d));
       bar(i, means(i), 'FaceColor', groups{i, 4}, 'EdgeColor', 'k', 'FaceAlpha', 0.8);
       text(i, means(i)+sems(i)+0.05*max(means), sprintf('%.2f', means(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'k');
   end
   errorbar(1:n_groups, means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
   xticks(1:n_groups); xticklabels(groups(:, 3)); xtickangle(15);
   y_start = max(means + sems) * 1.3; step = y_start * 0.15;
   for i = 2:n_groups
       tbl1 = groups{1, 1}; tbl2 = groups{i, 1};
       if ~isempty(tbl1) && ~isempty(tbl2)
           [~, p] = ttest2(tbl1.MedianLat, tbl2.MedianLat); if strcmp(mode, 'count'), [~, p] = ttest2(tbl1.NumResponses, tbl2.NumResponses); end
           h = y_start + (i-2)*step; plot([1, 1, i, i], [h-step*0.2, h, h, h-step*0.2], 'k-');
           text((1+i)/2, h + step*0.5, sprintf('p=%.1e', p), 'HorizontalAlignment', 'center', 'FontSize', 9);
       end
   end
   ylabel(lbl); title(figTitle); set(gca, 'FontSize', 11, 'Box', 'off'); hold off;
end

function plot_scatter_zoomed(groups, figTitle)
   figure('Name', figTitle, 'Color', 'w', 'Position', [50, 50, 1200, 400]);
   count = size(groups, 1); t = tiledlayout(1, count, 'Padding', 'compact', 'TileSpacing', 'compact');
   for i = 1:count
       nexttile; tbl = groups{i, 1}; col = groups{i, 4};
       if isempty(tbl), title([groups{i,3} ' (No Data)']); continue; end
       sel = 1 ./ (tbl.NumResponses ./ tbl.num_stimuli); lat = tbl.MedianLat;
       scatter(lat, sel, 15, 'filled', 'MarkerFaceColor', col, 'MarkerFaceAlpha', 0.7); hold on;
       l = lsline; l.Color = 'k'; l.LineWidth = 1;
       [rho_s, p_s] = corr(lat, sel, 'Type', 'Spearman'); [rho_p, p_p] = corr(lat, sel, 'Type', 'Pearson');
       stat_str = {sprintf('Spearman: \\rho=%.2f, p=%.2g', rho_s, p_s); sprintf('Pearson: \\rho=%.2f, p=%.2g', rho_p, p_p)};
       text(0.05, 0.95, stat_str, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9);
       title([groups{i, 3} ' (N=' num2str(height(tbl)) ')']); xlabel('latency [ms]'); if i==1, ylabel('selectivity'); end
       xlim([0 1000]); ylim([0 150]); box off; set(gca, 'TickDir', 'out'); hold off;
   end
   title(t, figTitle, 'FontSize', 14);
end

function plot_cdfs(groups, figTitle)
   figure('Name', figTitle, 'Color', 'w', 'Position', [100, 100, 500, 800]);
   t = tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
   nexttile; hold on; for i=1:size(groups,1), d=groups{i,1}.MedianLat; if ~isempty(d), h=cdfplot(d); set(h,'Color',groups{i,4},'LineWidth',2); end; end
   title('A. Latency'); xlabel('ms'); grid on; xlim([0 1000]); legend(groups(:,3),'Location','best'); hold off;
   nexttile; hold on; for i=1:size(groups,1), d=groups{i,1}.NumResponses; if ~isempty(d), h=cdfplot(d); set(h,'Color',groups{i,4},'LineWidth',2); end; end
   title('B. Responses/Unit'); xlabel('Count'); grid on; xlim([0 15]); hold off;
   nexttile; hold on; for i=1:size(groups,1), d=groups{i,1}.BaseFR; if ~isempty(d), h=cdfplot(d); set(h,'Color',groups{i,4},'LineWidth',2); end; end
   title('C. Firing Rate'); xlabel('Hz'); grid on; xlim([0 20]); hold off;
   delete(findall(gcf, 'Type', 'text', 'String', 'Empirical CDF'));
   title(t, figTitle, 'FontSize', 16);
end

function plot_combined_hist_bar(groups, figTitle, mode)
   figure('Name', figTitle, 'Color', 'w', 'Position', [50, 50, 1200, 800]);
   num_subregions = size(groups, 1);
   exclude_first_from_bars = contains(groups{1,3}, 'General', 'IgnoreCase', true);
   t = tiledlayout(2, num_subregions, 'Padding', 'compact', 'TileSpacing', 'compact');
   for i = 1:num_subregions
       nexttile(i); col = groups{i, 4}; name = groups{i, 3}; tbl = groups{i, 1};
       if isempty(tbl), title([name ' (No Data)']); continue; end
       if strcmp(mode, 'latency'), data = tbl.MedianLat; lbl='latency [ms]'; xlim_val = 1000; bin_args = {'BinWidth', 50}; fmt_str = '%.0f ms';
       elseif strcmp(mode, 'count'), data = tbl.NumResponses; lbl='responses per unit'; xlim_val = 15; bin_args = {'BinMethod', 'integers'}; fmt_str = '%.2f';
       elseif strcmp(mode, 'proportion'), data = tbl.NumResponses ./ tbl.num_stimuli; lbl='proportion of responses'; xlim_val = 1.0; bin_args = {'BinWidth', 0.05}; fmt_str = '%.3f'; end
       if isempty(data), title([name ' (No Data)']); continue; end
       h = histogram(data, bin_args{:}, 'FaceColor', col, 'FaceAlpha', 0.7, 'EdgeColor', 'k');
       maxY = max(h.Values); if maxY > 0, ylim([0, maxY * 1.2]); end
       med_val = median(data); xline(med_val, '--k', 'LineWidth', 1.5);
       title(sprintf('%s (N=%d)', name, height(tbl)), 'FontSize', 12); xlabel(lbl); if i==1, ylabel('number of units'); end
       xlim([0 xlim_val]); box on; set(gca, 'TickDir', 'out');
   end
   nexttile(num_subregions + 1, [1 num_subregions]); hold on;
   if exclude_first_from_bars, start_idx = 2; plot_indices = 2:num_subregions; else, start_idx = 1; plot_indices = 1:num_subregions; end
   num_bars = length(plot_indices);
   means = zeros(1, num_bars); sems = zeros(1, num_bars); all_data_for_ttest = cell(1, num_bars); bar_colors = cell(1, num_bars); bar_labels = cell(1, num_bars);
   for k = 1:num_bars
       real_idx = plot_indices(k); tbl = groups{real_idx, 1};
       if isempty(tbl), continue; end
       if strcmp(mode, 'latency'), d = tbl.MedianLat; elseif strcmp(mode, 'count'), d = tbl.NumResponses; elseif strcmp(mode, 'proportion'), d = tbl.NumResponses ./ tbl.num_stimuli; end
       if isempty(d), continue; end
       means(k) = mean(d); sems(k) = std(d)/sqrt(length(d)); all_data_for_ttest{k} = d; bar_colors{k} = groups{real_idx,4}; bar_labels{k} = groups{real_idx,3};
       bar(k, means(k), 'FaceColor', bar_colors{k}, 'EdgeColor', 'k', 'FaceAlpha', 0.8);
       text(k, means(k) + sems(k) + (0.1 * max(means)), sprintf(fmt_str, means(k)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 11, 'Color', 'k', 'FontWeight', 'bold');
   end
   errorbar(1:num_bars, means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
   xticks(1:num_bars); xticklabels(bar_labels); xtickangle(0); ylabel(lbl);
   if strcmp(mode, 'latency'), ylim([0 max(means)*2.0]); else, ylim([0 max(means)*2.5]); end
   box off; set(gca, 'TickDir', 'out', 'YGrid', 'on', 'GridAlpha', 0.15);
   yl = ylim; y_start_pval = yl(2) * 0.7; pval_step = yl(2) * 0.08;    
   if num_bars > 1 && ~isempty(all_data_for_ttest{1})
       for k = 2:num_bars
           if ~isempty(all_data_for_ttest{k})
               [~, p] = ttest2(all_data_for_ttest{1}, all_data_for_ttest{k});
               h_line = y_start_pval + (k-2)*pval_step;
               plot([1, 1, k, k], [h_line-pval_step*0.2, h_line, h_line, h_line-pval_step*0.2], 'k-');
               text((1+k)/2, h_line + pval_step*0.3, sprintf('p=%.1e', p), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
           end
       end
   end
   hold off; title(t, figTitle, 'FontSize', 16, 'FontWeight', 'bold');
end