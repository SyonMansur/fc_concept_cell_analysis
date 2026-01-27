function refinedTable = detectResponseProbabilistic(rawTable)
% parameters for state detection
bin_size = 20; % ms
min_consecutive = 3; % persistence requirement
analysis_window = [-500, 1000];

fprintf('refining latencies via hmm...\n');
numRows = height(rawTable);
rawTable.hmmResponseOnset = nan(numRows, 1);

% iterate through each response
for i = 1:numRows
    % check for spikes column added in extraction
    if ~ismember('spikes', rawTable.Properties.VariableNames)
        continue; 
    end
    
    spikeTimes = rawTable.spikes{i};
    allSpikes = vertcat(spikeTimes{:});
    
    % skip if silent
    if isempty(allSpikes)
        continue; 
    end

    % 1. bin the spikes
    time_bins = analysis_window(1):bin_size:analysis_window(2);
    counts = histcounts(allSpikes, time_bins);
    
    % 2. setup poisson-hmm states
    % state one: baseline, state two: response
    baseline_lambda = mean(counts(time_bins(1:end-1) < 0)) + 1e-6;
    response_lambda = mean(counts(time_bins(1:end-1) >= 0 & time_bins(1:end-1) < 400)) + 1e-5;
    
    % ensure response lambda is actually higher than baseline
    if response_lambda <= baseline_lambda
        response_lambda = baseline_lambda + 0.1;
    end
    
    % transition matrix: high probability of staying in current state
    trans_guess = [0.95, 0.05;  % p(1->1), p(1->2)
                   0.10, 0.90]; % p(2->1), p(2->2)
    
    % emission matrix: poisson probabilities for each possible spike count
    max_count = max(counts);
    obs_space = 0:max_count;
    emis_guess = [poisspdf(obs_space, baseline_lambda); 
                  poisspdf(obs_space, response_lambda)];
    emis_guess(emis_guess == 0) = 1e-10; % avoid zero probabilities
    
    try
        % 3. viterbi decoding
        % add one to counts because viterbi expects 1-based observations
        states = hmmviterbi(counts + 1, trans_guess, emis_guess);
        
        % 4. find transition with persistence
        stim_idx = find(time_bins >= 0, 1, 'first');
        post_stim_states = states(stim_idx:end);
        
        % look for the required number of consecutive response states
        target = repmat(2, 1, min_consecutive);
        match_idx = strfind(post_stim_states, target);
        
        if ~isempty(match_idx)
            % first index of the sequence is our onset
            onset_bin = stim_idx + match_idx(1) - 1;
            rawTable.hmmResponseOnset(i) = time_bins(onset_bin);
        end
    catch ME
        % ignore convergence failures and keep moving
        continue;
    end
end

refinedTable = rawTable;
fprintf('hmm refinement complete\n');
end