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

f5_chan = "F5";
target_chs = ["Fpz","AFz","Fz","FCz"];

for k = 1:numel(edf_pre)
    subj_pre = getSubjectID(edf_pre(k).name);
    participants{end+1} = subj_pre;

    % Match post file
    post_idx = find(arrayfun(@(x) strcmp(getSubjectID(x.name), subj_pre), edf_post), 1, 'first');

    % Load pre
    [pre_sig, pre_labels, pre_fs_all, pre_dur] = loadEDF(fullfile(pre_dir, edf_pre(k).name));
    fprintf('Subject %s PRE duration: %.2f sec\n', subj_pre, pre_dur);

    % Load post
    post_sig = []; post_labels = {}; post_fs_all = [];
    if ~isempty(post_idx)
        [post_sig, post_labels, post_fs_all, post_dur] = loadEDF(fullfile(post_dir, edf_post(post_idx).name));
        fprintf('Subject %s POST duration: %.2f sec\n', subj_pre, post_dur);
    end

    subj_row = [];

    % --- Compute bandpowers per channel ---
    for ch = 1:numel(channels)
        chan_name = channels(ch);

        % PRE
        [sig_pre, fs_pre] = extractChannel(pre_sig, pre_labels, pre_fs_all, chan_name);
        if ~isempty(sig_pre)
            bpp = computeBandPowersWelch(sig_pre, fs_pre, bands);
        else
            bpp = nan(1,4);
        end

        % POST
        if ~isempty(post_sig)
            [sig_post, fs_post] = extractChannel(post_sig, post_labels, post_fs_all, chan_name);
            bpp_post = computeBandPowersWelch(sig_post, fs_post, bands);
        else
            sig_post = []; fs_post = [];
            bpp_post = nan(1,4);
        end

        subj_row = [subj_row, bpp, bpp_post];
    end

    % --- PAC & PLV pre ---
    [preF5, fsF5] = extractChannel(pre_sig, pre_labels, pre_fs_all, f5_chan);
    pac_pre_theta_gamma = nan(1,numel(target_chs));
    pac_pre_theta_theta = nan(1,numel(target_chs));
    plv_pre_theta_theta = nan(1,numel(target_chs));

    for t = 1:numel(target_chs)
        tgt = target_chs(t);
        [preTgt, fsTgt] = extractChannel(pre_sig, pre_labels, pre_fs_all, tgt);
        if ~isempty(preF5) && ~isempty(preTgt)
            if fsF5 ~= fsTgt, preTgt = resample(preTgt, fsF5, fsTgt); fsTgt = fsF5; end
            pac_pre_theta_theta(t) = computeCrossPAC(preF5, fsF5, preTgt, fsTgt, [4 7], [4 7]);
            pac_pre_theta_gamma(t) = computeCrossPAC(preF5, fsF5, preTgt, fsTgt, [4 7], [31 45]);
            plv_pre_theta_theta(t) = computeCrossPLV(preF5, preTgt, fsF5, [4 7], [4 7]);
        end
    end
    subj_row = [subj_row, pac_pre_theta_gamma, pac_pre_theta_theta, plv_pre_theta_theta];

    % --- PAC & PLV post ---
    [postF5, fsF5_post] = extractChannel(post_sig, post_labels, post_fs_all, f5_chan);
    pac_post_theta_gamma = nan(1,numel(target_chs));
    pac_post_theta_theta = nan(1,numel(target_chs));
    plv_post_theta_theta = nan(1,numel(target_chs));

    for t = 1:numel(target_chs)
        tgt = target_chs(t);
        [postTgt, fsTgt_post] = extractChannel(post_sig, post_labels, post_fs_all, tgt);
        if ~isempty(postF5) && ~isempty(postTgt)
            if fsF5_post ~= fsTgt_post, postTgt = resample(postTgt, fsF5_post, fsTgt_post); fsTgt_post = fsF5_post; end
            pac_post_theta_theta(t) = computeCrossPAC(postF5, fsF5_post, postTgt, fsTgt_post, [4 7], [4 7]);
            pac_post_theta_gamma(t) = computeCrossPAC(postF5, fsF5_post, postTgt, fsTgt_post, [4 7], [31 45]);
            plv_post_theta_theta(t) = computeCrossPLV(postF5, postTgt, fsF5_post, [4 7], [4 7]);
        end
    end
    subj_row = [subj_row, pac_post_theta_gamma, pac_post_theta_theta, plv_post_theta_theta];

    results = [results; subj_row];
end

%% --- Build Table and Save ---
colnames = {};
for ch = channels
    colnames = [colnames, ...
        "Pre_"+ch+"_theta", "Pre_"+ch+"_alpha", "Pre_"+ch+"_beta", "Pre_"+ch+"_gamma", ...
        "Post_"+ch+"_theta", "Post_"+ch+"_alpha", "Post_"+ch+"_beta", "Post_"+ch+"_gamma"];
end

% PAC & PLV columns
for prefix = ["pre", "post"]
    for t = 1:numel(target_chs)
        chname = target_chs(t);
        colnames = [colnames, ...
            "PAC_"+prefix+"_thetaGamma_F5_"+chname, ...
            "PAC_"+prefix+"_thetaTheta_F5_"+chname, ...
            "PLV_"+prefix+"_thetaTheta_F5_"+chname];
    end
end

T = array2table(results, 'VariableNames', colnames);
T.Participant = participants';
T = movevars(T, 'Participant', 'before', 1);

writetable(T, 'EEG_features_crosscondition_prepost.csv');
fprintf('Saved EEG_features_crosscondition_prepost.csv\n');

%% ======================= HELPERS =======================

function subj = getSubjectID(fname)
    [~,name,~] = fileparts(fname);
    tokens = regexp(name, '_(p\d+[a-zA-Z]+)', 'tokens');
    if ~isempty(tokens), subj = tokens{1}{1}; else, subj = name; end
end

function [data, labels, fs_all, duration_sec] = loadEDF(filename)
    hdr = sopen(filename, 'r');
    sig = sread(hdr);
    labels = cellstr(hdr.Label); labels = strtrim(labels);
    fs_all = hdr.SampleRate;
    duration_sec = hdr.NRec * hdr.Dur;
    sclose(hdr);
    if size(sig,1) < size(sig,2), sig = sig'; end
    data = sig;
end

function [sig_col, fs_chan] = extractChannel(sig, labels, fs_all, target)
    if isempty(sig) || isempty(labels), sig_col = []; fs_chan = []; return; end
    idxs = find(strcmpi(labels, target));
    if isempty(idxs), sig_col = []; fs_chan = []; return; end
    idx = idxs(1);
    sig_col = double(sig(:, idx)); sig_col = sig_col(:);
    if isscalar(fs_all), fs_chan = fs_all; else, fs_chan = fs_all(idx); end
    sig_col = sig_col - mean(sig_col, 'omitnan');
end

function bp = computeBandPowersWelch(sig, fs, bands)
    if isempty(sig) || isempty(fs), bp = nan(1,numel(fieldnames(bands))); return; end
    winlen = min(4*fs, numel(sig)); if winlen < fs, winlen = numel(sig); end
    [pxx,f] = pwelch(sig, hamming(winlen), [], [], fs);
    fnames = fieldnames(bands); bp = nan(1,numel(fnames));
    for i = 1:numel(fnames)
        band = bands.(fnames{i});
        idx = f >= band(1) & f <= band(2);
        if nnz(idx) >= 2, bp(i) = trapz(f(idx), pxx(idx)); else, bp(i) = NaN; end
    end
end

function pac_val = computeCrossPAC(sigPre, fsPre, sigPost, fsPost, phaseBand, ampBand)
    if isempty(sigPre) || isempty(sigPost) || isempty(fsPre) || isempty(fsPost), pac_val = NaN; return; end
    if fsPre ~= fsPost, sigPost = resample(sigPost, fsPre, fsPost); fsPost = fsPre; end
    bp_phase = bandpass(sigPre, phaseBand, fsPre); phase = angle(hilbert(bp_phase));
    bp_amp = bandpass(sigPost, ampBand, fsPre); amp_env = abs(hilbert(bp_amp));
    nbins = 18; edges = linspace(-pi, pi, nbins+1);
    [~,~,bin] = histcounts(phase, edges);
    if all(bin==0), pac_val = NaN; return; end
    meanAmp = accumarray(bin(bin>0), amp_env(bin>0), [nbins 1], @mean, 0);
    if sum(meanAmp)==0, pac_val=NaN; return; end
    p = meanAmp / sum(meanAmp); H = -nansum(p .* log(p + eps)); Hmax = log(nbins);
    pac_val = (Hmax - H) / Hmax;
end

function plv = computeCrossPLV(sigPre, sigPost, fs, bandPre, bandPost)
    if isempty(sigPre) || isempty(sigPost) || isempty(fs), plv = NaN; return; end
    xpre = bandpass(sigPre, bandPre, fs); xpost = bandpass(sigPost, bandPost, fs);
    phi1 = angle(hilbert(xpre)); phi2 = angle(hilbert(xpost));
    n = min(numel(phi1), numel(phi2)); if n==0, plv=NaN; return; end
    plv = abs(mean(exp(1i*(phi1(1:n) - phi2(1:n)))));
end