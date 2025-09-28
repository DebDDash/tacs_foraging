
data = readtable("EEG_features_crosscondition_prepost.csv");

% Electrode pairs of interest
pairs = {'F5_Fpz','F5_FCz','F5_AFz','F5_Fz'};

% Initialize results
pre_means = zeros(1, numel(pairs));
post_means = zeros(1, numel(pairs));

% Loop over electrode pairs
for i = 1:numel(pairs)
    pre_col = sprintf("PAC_pre_thetaGamma_%s", pairs{i});
    post_col = sprintf("PAC_post_thetaGamma_%s", pairs{i});
    
    % Compute group average
    pre_means(i) = mean(data.(pre_col), 'omitnan');
    post_means(i) = mean(data.(post_col), 'omitnan');
end

% Prepare matrix for heatmap (2 x 4: Pre vs Post × Electrode pairs)
heatmap_matrix = [pre_means; post_means];

% Create heatmap
figure;
h = heatmap(pairs, {'Pre','Post'}, heatmap_matrix, ...
            'Colormap', parula, 'ColorbarVisible','on');
h.Title = 'Group-Averaged PAC Theta–Gamma Coupling';
h.XLabel = 'Electrode Pair';
h.YLabel = 'Condition';


