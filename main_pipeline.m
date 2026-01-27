%% fried lab data pipeline
% Coordinates the flow from raw clusterCharacteristics.mat files to final figures
% UPDATED 12/11/2025


% paths!

% directory of all raw data should go here!
rawDir = '';

% directory of all localization data should go here!
locDir = '';

% final output excel sheet here
output_xlsx = '';

%% 1. extraction
% loads clusterCharacteristics.mat files and runs initial zeta-test
% script: extractSelectiveResponses.m
fprintf('step one: extracting selective responses\n');
raw_table = extractSelectiveResponses(raw_data_path);

%% 2. probabilistic refinement (HMM)
% refines response onsets using a two-state poisson hidden markov model
% script: detectResponseProbabilistic.m
fprintf('step two: refining latencies via hmm\n');
% this uses the viterbi algorithm to find state transitions
raw_table = detectResponseProbabilistic(raw_table);

%% 3. regional mapping
% adds cleaned anatomical labels and broad areas (MTL/FC)
% script: mapBrainRegions.m
fprintf('step three: mapping brain regions\n');
mapped_table = mapBrainRegions(raw_table);

%% 4. channel filtering
% selects the "best" cluster per channel based on highest mean score
% script: filterBestClusters.m
fprintf('step four: filtering for best clusters\n');
final_table = filterBestClusters(mapped_table);

%% 5. anatomical verification
% cross-references with excel localization sheets for "trueLocation"
% script: verifyLocalizations.m
fprintf('step five: verifying localizations\n');
final_table = verifyLocalizations(final_table, localization_dir);

%% 6. visualization
% generates all histograms, bar plots, and scatters
% script: generatePipelineFigures.m
fprintf('step six: generating figures\n');
generatePipelineFigures(final_table);

%% 7. export
% saves the processed table to excel
fprintf('saving results to %s\n', output_xlsx);
writetable(final_table, output_xlsx);

fprintf('pipeline complete\n');