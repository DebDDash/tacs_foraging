%% --------------------- Setup ---------------------
% ---- BIOSIG path ----
biosig_dir = '/Users/debarpita/Desktop/arjun/biosig4octmat-3.8.5/biosig';
addpath(genpath(biosig_dir));

% --- Directories ---
pre_dir  = '/Users/debarpita/Desktop/arjun/pre';
post_dir = '/Users/debarpita/Desktop/arjun/post';

edf_pre  = dir(fullfile(pre_dir, '*.edf'));
edf_post = dir(fullfile(post_dir, '*.edf'));

% --- Parameters ---
channels = ["Fz","FCz","F5","Fpz","AFz"];  
phase_freqs = 2:2:20;     % low-frequency phase (Hz)
amp_freqs   = 30:10:150;  % high-frequency amplitude (Hz)
nBins = 18;               % number of phase bins

%% --------------------- Helper Functions ---------------------
function subjID = getSubjectID(filename)
    % Extract subject ID assuming pattern 'sub-XX_...'
    tokens = regexp(filename,'sub-(\d+)','tokens');
    if ~isempty(tokens)
        subjID = tokens{1}{1};
    else
        subjID = filename;
    end
end

function [sig, labels, fs, duration] = loadEDF(file)
    % Load EDF file using BIOSIG
    hdr = sopen(file);
    fs = hdr.SampleRate(1);
    nSamples = hdr.NS;
    duration = nSamples/fs;
    labels = hdr.Label;
    sig = sread(hdr)';  % transpose to channels x time
    sclose(hdr);
end

%% --------------------- Main Loop ---------------------
participants = {};
comod_all = [];  % store comodulograms for all subjects

for k = 1:numel(edf_pre)
    % --- Subject ID ---
    subj_pre = getSubjectID(edf_pre(k).name);
    participants{end+1} = subj_pre;

    % --- Load PRE data ---
    [pre_sig, pre_labels, fs, pre_dur] = loadEDF(fullfile(pre_dir, edf_pre(k).name));
    fprintf('Subject %s PRE duration: %.2f sec\n', subj_pre, pre_dur);

    % Select channels of interest
    chan_idx = find(ismember(pre_labels, channels));
    data = pre_sig(chan_idx,:);

    nPhase = numel(phase_freqs);
    nAmp = numel(amp_freqs);
    nChan = size(data,1);

    comod_subj = zeros(nAmp, nPhase, nChan);

    % --- Compute PAC ---
    for ch = 1:nChan
        sig_ch = data(ch,:);

        for ip = 1:nPhase
            pf = phase_freqs(ip);
            phase_filt = bandpass(sig_ch, [pf-1 pf+1], fs);
            phase_sig = angle(hilbert(phase_filt));

            for ia = 1:nAmp
                af = amp_freqs(ia);
                amp_filt = bandpass(sig_ch, [af-5 af+5], fs);
                amp_sig = abs(hilbert(amp_filt));

                % Modulation Index (Tort et al., 2010)
                edges = linspace(-pi, pi, nBins+1);
                [~,~,bin] = histcounts(phase_sig, edges);
                bin(bin==0) = 1;  % assign out-of-range samples to first bin
                meanAmp = accumarray(bin(:), amp_sig(:), [nBins 1], @mean, 0);
                meanAmp = meanAmp / sum(meanAmp);
                H = -sum(meanAmp .* log(meanAmp+eps));
                MI = (log(nBins) - H)/log(nBins);

                comod_subj(ia,ip,ch) = MI;
            end
        end
    end

    % Average across channels
    comod_all(:,:,k) = mean(comod_subj,3,'omitnan');
end

%% --------------------- Group-level PAC ---------------------
comod_group = mean(comod_all,3,'omitnan');

%% --------------------- Plotting ---------------------
figure;
imagesc(phase_freqs, amp_freqs, comod_group);
set(gca,'YDir','normal');
xlabel('Phase Frequency (Hz)');
ylabel('Amplitude Frequency (Hz)');
title('Group-level PAC Comodulogram');
colormap(jet);
colorbar;

for k = 1:numel(edf_post)
    % --- Subject ID ---
    subj_post = getSubjectID(edf_post(k).name);
    participants{end+1} = subj_post;

    % --- Load PRE data ---
    [post_sig, post_labels, fs, post_dur] = loadEDF(fullfile(post_dir, edf_post(k).name));
    fprintf('Subject %s PRE duration: %.2f sec\n', subj_pre, post_dur);

    % Select channels of interest
    chan_idx = find(ismember(post_labels, channels));
    data = post_sig(chan_idx,:);

    nPhase = numel(phase_freqs);
    nAmp = numel(amp_freqs);
    nChan = size(data,1);

    comod_subj = zeros(nAmp, nPhase, nChan);

    % --- Compute PAC ---
    for ch = 1:nChan
        sig_ch = data(ch,:);

        for ip = 1:nPhase
            pf = phase_freqs(ip);
            phase_filt = bandpass(sig_ch, [pf-1 pf+1], fs);
            phase_sig = angle(hilbert(phase_filt));

            for ia = 1:nAmp
                af = amp_freqs(ia);
                amp_filt = bandpass(sig_ch, [af-5 af+5], fs);
                amp_sig = abs(hilbert(amp_filt));

                % Modulation Index (Tort et al., 2010)
                edges = linspace(-pi, pi, nBins+1);
                [~,~,bin] = histcounts(phase_sig, edges);
                bin(bin==0) = 1;  % assign out-of-range samples to first bin
                meanAmp = accumarray(bin(:), amp_sig(:), [nBins 1], @mean, 0);
                meanAmp = meanAmp / sum(meanAmp);
                H = -sum(meanAmp .* log(meanAmp+eps));
                MI = (log(nBins) - H)/log(nBins);

                comod_subj(ia,ip,ch) = MI;
            end
        end
    end

    % Average across channels
    comod_all(:,:,k) = mean(comod_subj,3,'omitnan');
end

%% --------------------- Group-level PAC ---------------------
comod_group = mean(comod_all,3,'omitnan');

%% --------------------- Plotting ---------------------
figure;
imagesc(phase_freqs, amp_freqs, comod_group);
set(gca,'YDir','normal');
xlabel('Phase Frequency (Hz)');
ylabel('Amplitude Frequency (Hz)');
title('Group-level PAC Comodulogram');
colormap(jet);
colorbar;