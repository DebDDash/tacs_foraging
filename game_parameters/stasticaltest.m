%% Statistical analysis for EEG metrics
clear; clc;

% Load data
data = readtable('EEG_metrics_all_participants.csv');

% Metrics to test (must match your CSV column names)
metrics = {'ThetaBeta', 'AlphaBeta', 'PAC', 'Alpha', 'Beta', 'Gamma'};

% Initialize results table
Results = table('Size', [length(metrics), 7], ...
    'VariableTypes', {'string','double','double','double','double','double','double'}, ...
    'VariableNames', {'Metric','PreMean','PostMean','tStat','pValue','CohensD','Wilcoxon_p'});

for i = 1:length(metrics)
    m = metrics{i};
    pre_col = sprintf('Pre_%s', m);
    post_col = sprintf('Post_%s', m);

    pre_vals = data.(pre_col);
    post_vals = data.(post_col);

    % Paired t-test
    [~, p, ~, stats] = ttest(pre_vals, post_vals);
    tval = stats.tstat;

    % Effect size (Cohen's d for paired)
    diff_vals = post_vals - pre_vals;
    d = mean(diff_vals) / std(diff_vals);

    % Nonparametric alternative (Wilcoxon signed-rank)
    p_wilcoxon = signrank(pre_vals, post_vals);

    % Mean values
    pre_mean = mean(pre_vals);
    post_mean = mean(post_vals);

    % Store
    Results(i,:) = {m, pre_mean, post_mean, tval, p, d, p_wilcoxon};
end

% Save results to CSV
writetable(Results, 'EEG_pre_post_stats.csv');

% Display results
disp('✅ Statistical test results:')
disp(Results)