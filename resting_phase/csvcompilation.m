% ---- BIOSIG path ----
biosig_dir = '/Users/debarpita/Desktop/arjun/biosig4octmat-3.8.5/biosig';
addpath(genpath(biosig_dir));

%% --- Directories ---
pre_dir  = '/Users/debarpita/Desktop/arjun/pre';
post_dir = '/Users/debarpita/Desktop/arjun/post';

edf_pre  = dir(fullfile(pre_dir, '*.edf'));
edf_post = dir(fullfile(post_dir, '*.edf'));

%% --- Parameters ---
channels = ["Fz","FCz","F5","Fpz","AFz"];  
bands = struct('theta',[4 7], 'alpha',[8 12], 'beta',[13 30], 'gamma',[31 45]);

%% --- Main loop ---
participants = {};
results = [];

for k = 1:numel(edf_pre)
    subj_pre = getSubjectID(edf_pre(k).name);
    participants{end+1} = subj_pre;

    % Match post file (first match)
    post_idx = find(arrayfun(@(x) strcmp(getSubjectID(x.name), subj_pre), edf_post), 1, 'first');

    % Load pre
    [pre_sig, pre_labels, pre_fs_all, pre_dur] = loadEDF(fullfile(pre_dir, edf_pre(k).name));
    fprintf('Subject %s PRE duration: %.2f sec\n', subj_pre, pre_dur);

    % Load post (if exists)
    post_sig = []; post_labels = {}; post_fs_all = [];
    if ~isempty(post_idx)
        [post_sig, post_labels, post_fs_all, post_dur] = loadEDF(fullfile(post_dir, edf_post(post_idx).name));
        fprintf('Subject %s POST duration: %.2f sec\n', subj_pre, post_dur);
    end

    subj_row = [];

    % --- For every channel: compute pre & post bandpowers (4), and
    % --- cross-condition PACs: pre-theta phase vs post-theta amp & post-gamma amp
    for ch = 1:numel(channels)
        chan_name = channels(ch);

        % --- PRE: extract channel
        [sig_pre, fs_pre] = extractChannel(pre_sig, pre_labels, pre_fs_all, chan_name);
        if ~isempty(sig_pre)
            bpp = computeBandPowersWelch(sig_pre, fs_pre, bands);
        else
            bpp = nan(1,4);
        end

        % --- POST: extract channel
        if ~isempty(post_sig)
            [sig_post, fs_post] = extractChannel(post_sig, post_labels, post_fs_all, chan_name);
        else
            sig_post = []; fs_post = [];
        end

        % Append bandpowers: pre then post
        subj_row = [subj_row, bpp];
        if ~isempty(sig_post)
            bpp_post = computeBandPowersWelch(sig_post, fs_post, bands);
            subj_row = [subj_row, bpp_post];
        else
            subj_row = [subj_row, nan(1,4)];
        end

        % --- Cross-condition PAC for this channel:
        % phase from PRE theta; amplitude from POST theta (theta-theta)
        % and amplitude from POST gamma (theta-gamma)
        if ~isempty(sig_pre) && ~isempty(sig_post)
            % Ensure same fs: resample post to pre fs if necessary
            if isempty(fs_pre), fs_pre = fs_post; end
            if isempty(fs_post), fs_post = fs_pre; end
            if fs_post ~= fs_pre
                sig_post_rs = resample(sig_post, fs_pre, fs_post);
                fs_post_used = fs_pre;
            else
                sig_post_rs = sig_post;
                fs_post_used = fs_post;
            end

            pac_tt = computeCrossPAC(sig_pre, fs_pre, sig_post_rs, fs_post_used, [4 7], [4 7]);
            pac_tg = computeCrossPAC(sig_pre, fs_pre, sig_post_rs, fs_post_used, [4 7], [31 45]);
            subj_row = [subj_row, pac_tg, pac_tt];  % keep order theta-gamma then theta-theta
        else
            subj_row = [subj_row, nan(1,2)];
        end
    end

    % --- Cross-condition PLVs:
    % Use Pre F5 phase vs Post (Fpz, AFz, Fz, FCz) for:
    % 1) pre-theta vs post-theta PLV
    % 2) pre-theta vs post-gamma PLV
    plv_list = nan(1,8); % 4 pairs x 2 metrics
    target_chs = ["Fpz","AFz","Fz","FCz"];
    % Extract pre F5
    [preF5, preF5_fs] = extractChannel(pre_sig, pre_labels, pre_fs_all, "F5");
    if isempty(preF5)
        % cannot compute cross PLVs -> remain NaNs
        subj_row = [subj_row, plv_list];
    else
        % For each target channel compute both PLV types
        col = 1;
        for t = 1:numel(target_chs)
            tgt = target_chs(t);
            if ~isempty(post_sig)
                [postSigTgt, post_fs] = extractChannel(post_sig, post_labels, post_fs_all, tgt);
            else
                postSigTgt = []; post_fs = [];
            end

            if isempty(postSigTgt)
                plv_theta_theta = NaN;
                plv_theta_gamma = NaN;
            else
                % resample post target to preF5 fs if needed
                if isempty(preF5_fs), preF5_fs = post_fs; end
                if isempty(post_fs), post_fs = preF5_fs; end
                if preF5_fs ~= post_fs
                    post_rs = resample(postSigTgt, preF5_fs, post_fs);
                else
                    post_rs = postSigTgt;
                end

                % PLV pre-theta vs post-theta
                plv_theta_theta = computeCrossPLV(preF5, post_rs, preF5_fs, [4 7], [4 7]);

                % PLV pre-theta vs post-gamma (phase-of-gamma)
                plv_theta_gamma = computeCrossPLV(preF5, post_rs, preF5_fs, [4 7], [31 45]);
            end

            plv_list(col) = plv_theta_theta;   col = col + 1;
            plv_list(col) = plv_theta_gamma;   col = col + 1;
        end
        subj_row = [subj_row, plv_list];
    end

    results = [results; subj_row];
end

%% --- Build Table and Save ---
colnames = {};
for ch = channels
    % For each channel: Pre bandpowers (4), Post bandpowers (4), PAC_tg, PAC_tt
    colnames = [colnames, ...
        "Pre_"+ch+"_theta", "Pre_"+ch+"_alpha", "Pre_"+ch+"_beta", "Pre_"+ch+"_gamma", ...
        "Post_"+ch+"_theta","Post_"+ch+"_alpha","Post_"+ch+"_beta","Post_"+ch+"_gamma", ...
        "CrossPAC_"+ch+"_thetaGamma", ...
        "CrossPAC_"+ch+"_thetaTheta"];
end

% PLV columns: for each target channel, two metrics (theta-theta, theta-gamma)
target_chs = ["Fpz","AFz","Fz","FCz"];
for t = 1:numel(target_chs)
    chname = target_chs(t);
    colnames = [colnames, ...
        "CrossPLV_PreF5_Post"+chname+"_theta_theta", ...
        "CrossPLV_PreF5_Post"+chname+"_theta_gamma"];
end

% final table
T = array2table(results, 'VariableNames', colnames);
T.Participant = participants';
T = movevars(T, 'Participant', 'before', 1);

writetable(T, 'EEG_features_crosscondition_prepost.csv');
fprintf('Saved EEG_features_crosscondition_prepost.csv\n');

%% ======================= HELPERS =======================

function subj = getSubjectID(fname)
    [~,name,~] = fileparts(fname);
    tokens = regexp(name, '_(p\d+[a-zA-Z]+)', 'tokens');
    if ~isempty(tokens)
        subj = tokens{1}{1};
    else
        subj = name;
    end
end

function [data, labels, fs_all, duration_sec] = loadEDF(filename)
    hdr = sopen(filename, 'r');
    sig = sread(hdr);
    labels = cellstr(hdr.Label); labels = strtrim(labels);
    fs_all = hdr.SampleRate;
    duration_sec = hdr.NRec * hdr.Dur;
    sclose(hdr);
    if size(sig,1) < size(sig,2)
        sig = sig';
    end
    data = sig;
end

function [sig_col, fs_chan] = extractChannel(sig, labels, fs_all, target)
    % safe extractor: returns column vector (samples x 1) and channel fs
    if isempty(sig) || isempty(labels)
        sig_col = []; fs_chan = [];
        return;
    end
    idxs = find(strcmpi(labels, target));
    if isempty(idxs)
        sig_col = []; fs_chan = [];
        return;
    end
    idx = idxs(1);
    sig_col = double(sig(:, idx));
    if isscalar(fs_all)
        fs_chan = fs_all;
    else
        fs_chan = fs_all(idx);
    end
    sig_col = sig_col(:);
    sig_col = sig_col - mean(sig_col, 'omitnan');
end

function bp = computeBandPowersWelch(sig, fs, bands)
    if isempty(sig) || isempty(fs)
        bp = nan(1,numel(fieldnames(bands))); return;
    end
    winlen = min(4*fs, numel(sig));
    if winlen < fs, winlen = numel(sig); end
    [pxx,f] = pwelch(sig, hamming(winlen), [], [], fs);
    fnames = fieldnames(bands);
    bp = nan(1,numel(fnames));
    for i = 1:numel(fnames)
        band = bands.(fnames{i});
        idx = f >= band(1) & f <= band(2);
        if nnz(idx) >= 2
            bp(i) = trapz(f(idx), pxx(idx));
        else
            bp(i) = NaN;
        end
    end
end

function pac_val = computeCrossPAC(sigPre, fsPre, sigPost, fsPost, phaseBand, ampBand)
    % pre: phaseBand -> get phase
    % post: ampBand -> get amplitude envelope
    if isempty(sigPre) || isempty(sigPost) || isempty(fsPre) || isempty(fsPost)
        pac_val = NaN; return;
    end
    % Ensure signals are same fs (sigPost should already be resampled upstream if needed)
    if fsPre ~= fsPost
        sigPost = resample(sigPost, fsPre, fsPost);
        fsPost = fsPre;
    end

    % Bandpass and Hilbert
    bp_phase = bandpass(sigPre, phaseBand, fsPre);
    phase = angle(hilbert(bp_phase));

    bp_amp = bandpass(sigPost, ampBand, fsPre);
    amp_env = abs(hilbert(bp_amp));

    % compute Modulation Index (Tort)
    nbins = 18;
    edges = linspace(-pi, pi, nbins+1);
    [~,~,bin] = histcounts(phase, edges);
    if all(bin==0)
        pac_val = NaN; return;
    end
    meanAmp = accumarray(bin(bin>0), amp_env(bin>0), [nbins 1], @mean, 0);
    if sum(meanAmp) == 0
        pac_val = NaN; return;
    end
    p = meanAmp / sum(meanAmp);
    H = -nansum(p .* log(p + eps));
    Hmax = log(nbins);
    pac_val = (Hmax - H) / Hmax;
end

function plv = computeCrossPLV(sigPre, sigPost, fs, bandPre, bandPost)
    % Compute PLV between phase(sigPre in bandPre) and phase(sigPost in bandPost)
    if isempty(sigPre) || isempty(sigPost) || isempty(fs)
        plv = NaN; return;
    end
    % bandpass
    xpre = bandpass(sigPre, bandPre, fs);
    xpost = bandpass(sigPost, bandPost, fs);
    phi1 = angle(hilbert(xpre));
    phi2 = angle(hilbert(xpost));
    % ensure same length
    n = min(numel(phi1), numel(phi2));
    if n == 0
        plv = NaN; return;
    end
    phi1 = phi1(1:n);
    phi2 = phi2(1:n);
    plv = abs(mean(exp(1i*(phi1 - phi2))));
end