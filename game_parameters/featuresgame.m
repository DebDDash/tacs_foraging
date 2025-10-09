biosig_dir = '/Users/debarpita/Desktop/arjun/biosig4octmat-3.8.5/biosig';
eeglab_dir = '/Users/debarpita/Desktop/arjun/eeglab2025.0.0';
addpath(genpath(biosig_dir));
addpath(genpath(eeglab_dir));

%% Initialize EEGLAB
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Define participants
participants = {'p1', 'p2', 'p3', 'p4','p5','p6','p7','p8','p9','p10'};  % add all participants here
results = [];

% Define frequency bands (Hz)
theta_band = [4 7];
alpha_band = [8 13];
beta_band  = [14 30];
gamma_band = [31 45];

% Define PAC frequency bands (example: theta phase and gamma amplitude)
pac_phase_band = [4 7];
pac_amp_band   = [31 45];

function [theta_pow, alpha_pow, beta_pow, gamma_pow] = band_power(EEG, theta_band, alpha_band, beta_band, gamma_band)
    % Compute mean power in frequency bands using pwelch
    fs = EEG.srate;
    data = EEG.data;
    n_channels = size(data, 1);
    
    freq_res = [];
    for ch = 1:n_channels
        [pxx, f] = pwelch(data(ch,:), [], [], [], fs);
        freq_res(ch,:) = pxx';
    end
    
    % Mean over channels
    mean_pxx = mean(freq_res, 1);

    theta_pow = bandpower(mean_pxx, f, theta_band, 'psd');
    alpha_pow = bandpower(mean_pxx, f, alpha_band, 'psd');
    beta_pow  = bandpower(mean_pxx, f, beta_band, 'psd');
    gamma_pow = bandpower(mean_pxx, f, gamma_band, 'psd');
end

function pac_val = compute_PAC(EEG, phase_band, amp_band)
    data = mean(EEG.data, 1);  % average across channels
    fs = EEG.srate;

    % Bandpass filter for phase and amplitude
    phase_sig = eegfilt(data, fs, phase_band(1), phase_band(2));
    amp_sig   = eegfilt(data, fs, amp_band(1), amp_band(2));

    % Get phase and amplitude envelope
    phase = angle(hilbert(phase_sig));
    amp_env = abs(hilbert(amp_sig));

    % Bin the phase into 18 bins
    n_bins = 18;
    edges = linspace(-pi, pi, n_bins+1);
    mean_amp = zeros(1, n_bins);

    for b = 1:n_bins
        idx = phase >= edges(b) & phase < edges(b+1);
        mean_amp(b) = mean(amp_env(idx));
    end

    mean_amp = mean_amp / sum(mean_amp);
    H = -sum(mean_amp .* log(mean_amp + eps));
    Hmax = log(n_bins);
    pac_val = (Hmax - H) / Hmax;  % Modulation Index (0–1)
end


% Loop through each participant
for i = 1:length(participants)
    subj = participants{i};

    % File paths
    pre_file  = sprintf('%s_prestim_epochs.set', subj);
    post_file = sprintf('%s_poststim_epochs.set', subj);

    %% ----- PRESTIM -----
    EEG_pre = pop_loadset('filename', pre_file);
    EEG_pre = eeg_checkset(EEG_pre);
    [pre_theta, pre_alpha, pre_beta, pre_gamma] = band_power(EEG_pre, theta_band, alpha_band, beta_band, gamma_band);

    % Ratios
    pre_theta_beta_ratio = pre_theta / pre_beta;
    pre_alpha_beta_ratio = pre_alpha / pre_beta;

    % PAC (Phase-Amplitude Coupling)
    pre_pac = compute_PAC(EEG_pre, pac_phase_band, pac_amp_band);

    %% ----- POSTSTIM -----
    EEG_post = pop_loadset('filename', post_file);
    EEG_post = eeg_checkset(EEG_post);
    [post_theta, post_alpha, post_beta, post_gamma] = band_power(EEG_post, theta_band, alpha_band, beta_band, gamma_band);

    post_theta_beta_ratio = post_theta / post_beta;
    post_alpha_beta_ratio = post_alpha / post_beta;
    post_pac = compute_PAC(EEG_post, pac_phase_band, pac_amp_band);

    %% Save results for this subject
    results = [results; {subj, ...
        pre_theta_beta_ratio, post_theta_beta_ratio, ...
        pre_alpha_beta_ratio, post_alpha_beta_ratio, ...
        pre_pac, post_pac, ...
        pre_alpha, post_alpha, pre_beta, post_beta, pre_gamma, post_gamma}];
end

%% Convert to table and save as CSV
headers = {'Participant', ...
    'Pre_ThetaBeta', 'Post_ThetaBeta', ...
    'Pre_AlphaBeta', 'Post_AlphaBeta', ...
    'Pre_PAC', 'Post_PAC', ...
    'Pre_Alpha', 'Post_Alpha', 'Pre_Beta', 'Post_Beta', 'Pre_Gamma', 'Post_Gamma'};

T = cell2table(results, 'VariableNames', headers);
writetable(T, 'EEG_metrics_all_participants.csv');

disp('✅ Analysis complete! Results saved to EEG_metrics_all_participants.csv');