%% This code was repeated for alpha,beta and gamma range
%% --- Setup BIOSIG Path ---
biosig_dir = '/Users/debarpita/Desktop/arjun/biosig4octmat-3.8.5/biosig';
addpath(genpath(biosig_dir));

%% --- Directories for pre and post ---
pre_dir  = '/Users/debarpita/Desktop/arjun/pre';
post_dir = '/Users/debarpita/Desktop/arjun/post';

edf_pre  = dir(fullfile(pre_dir, '*.edf'));
edf_post = dir(fullfile(post_dir, '*.edf'));

%% --- Parameters ---
channels = ["Fpz","F5","Fcz","Fz"]; % add channels you want
epochLength = 2; % seconds per epoch

%% --- Helper: Extract subject ID from filename ---
%% --- Helper: Extract subject ID from filename ---
function subj = getSubjectID(fname)
    [~,name,~] = fileparts(fname);
    % keep only the "p4ashwin" part, ignore session/run suffixes
    tokens = regexp(name, '_(p\d+[a-zA-Z]+)', 'tokens');
    if ~isempty(tokens)
        subj = tokens{1}{1};
    else
        subj = name; % fallback if regex fails
    end
end

%% --- Helper: Power computation ( in this case we consider Low gamma band 30–45 Hz) ---
function powerVals = computePower(sig, fs, samplesPerEpoch)
    % Bandpass low gamma (30–45 Hz)
    bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
        'HalfPowerFrequency1',30,'HalfPowerFrequency2',45, ...
        'SampleRate',fs);
    sig_gamma = filtfilt(bpFilt, sig);
    % Per-epoch variance as power
    powerVals = var(buffer(sig_gamma, samplesPerEpoch, 0, 'nodelay'));
end

%% --- Process files function ---
function T = processFiles(edf_list, data_dir, channels, cond, epochLength)
    T = table();
    for i = 1:length(edf_list)
        fname = fullfile(data_dir, edf_list(i).name);
        hdr = sopen(fname); [data,hdr] = sread(hdr); sclose(hdr);
        fs = hdr.SampleRate;
        samplesPerEpoch = epochLength * fs;

        for c = channels
            ch_idx = find(strcmpi(hdr.Label, c),1);
            if isempty(ch_idx)
                warning(['Channel ',c,' not found in ',fname]);
                continue;
            end
            sig = data(:,ch_idx);
            powerVals = computePower(sig, fs, samplesPerEpoch);

            tmp = table();
            tmp.Subject    = repmat(string(getSubjectID(fname)), numel(powerVals),1);
            tmp.Channel    = repmat(string(c), numel(powerVals),1);
            tmp.Condition  = repmat(string(cond), numel(powerVals),1);
            tmp.ThetaPower = powerVals(:);
            T = [T; tmp];
        end
    end
end

%% --- Load all data ---
T_pre  = processFiles(edf_pre,  pre_dir,  channels, 'Pre', epochLength);
T_post = processFiles(edf_post, post_dir, channels, 'Post', epochLength);
T = [T_pre; T_post];

%% --- Combined boxplot ---
figure;
boxplot(T.ThetaPower, {T.Channel, T.Condition});
ylabel('Gamma Power (\muV^2)');
title('Low Gamma Power: Pre vs Post Across All Subjects');

%% --- Statistical Analysis ---
%% --- Statistical Analysis ---
comparisons = {
    {'Fz','Pre','Post','Fz Pre vs Post'}
    {'Fpz','Pre','Post','Fpz Pre vs Post'}
    {'F5','Pre','Fpz','Post','F5 Pre vs Fpz Post'}
    {'F5','Pre','Fcz','Post','F5 Pre vs Fcz Post'}
    {'Fcz','Pre','Post','Fcz Pre vs Post'}
    {'Fpz','Pre','Fcz','Post','Fpz Pre vs Fcz Post'}
};

Results = table();

for i = 1:length(comparisons)
    comp = comparisons{i};

    if numel(comp)==4
        % same channel pre vs post
        ch = comp{1}; condA = comp{2}; condB = comp{3}; label = comp{4};
        subjList = intersect(unique(T.Subject(T.Channel==ch & T.Condition==condA)), ...
                             unique(T.Subject(T.Channel==ch & T.Condition==condB)));
        x=[]; y=[];
        for s = subjList'
            xa = mean(T.ThetaPower(T.Subject==s & T.Channel==ch & T.Condition==condA));
            yb = mean(T.ThetaPower(T.Subject==s & T.Channel==ch & T.Condition==condB));
            x(end+1)=xa; y(end+1)=yb;
        end
    else
        % cross channel comparisons
        chA = comp{1}; condA = comp{2}; chB = comp{3}; condB = comp{4}; label = comp{5};
        subjList = intersect(unique(T.Subject(T.Channel==chA & T.Condition==condA)), ...
                             unique(T.Subject(T.Channel==chB & T.Condition==condB)));
        x=[]; y=[];
        for s = subjList'
            xa = mean(T.ThetaPower(T.Subject==s & T.Channel==chA & T.Condition==condA));
            yb = mean(T.ThetaPower(T.Subject==s & T.Channel==chB & T.Condition==condB));
            x(end+1)=xa; y(end+1)=yb;
        end
    end

    % ✅ Debug print goes HERE
    fprintf('%s: n = %d subjects\n', label, numel(x));
    disp(table(subjList, x', y', (y-x)'))

    if isempty(x) || isempty(y)
        continue;
    end

    % Wilcoxon signed-rank
    % Wilcoxon signed-rank
    try
        [p,~,stats] = signrank(x,y);
        zval = stats.zval;
    catch
    % fallback: paired t-test
    [~,p,~,stats] = ttest(x,y);
    zval = stats.tstat;  % use t instead of z
    end

    % Cohen's d
    d = (mean(y-x)) / std(y-x);

    % Effect size r (non-parametric)
    r = zval / sqrt(length(x));

    meanDiff = mean(y-x);
    pctChange = 100 * meanDiff ./ mean(x);

    newRow = table(string(label), p, zval, d, r, meanDiff, pctChange, ...
    'VariableNames', {'Comparison','pValue','zValue','CohensD','EffectSizeR','MeanDiff','PctChange'});
    Results = [Results; newRow];
end

disp(Results);
