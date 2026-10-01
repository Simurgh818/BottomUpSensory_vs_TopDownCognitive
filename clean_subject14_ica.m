%% =========================================================================
%% CLEAN SUBJECT 14 VIA ICA & VERIFY ARTIFACT REMOVAL
%% =========================================================================
clear; clc; close all;

if exist('I:\', 'dir')
    input_path = 'I:\My Drive\Data\New Data\EEG epoched\';
    base_output_path = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\Dr. Sederberg MaTRIX Lab\Research Paper\';
elseif exist('H:\', 'dir')
    input_path = 'H:\My Drive\Data\New Data\EEG epoched\';
    base_output_path = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\Dr. Sederberg MaTRIX Lab\Research Paper\';
elseif exist('G:\', 'dir')
    input_path = 'G:\My Drive\Data\New Data\EEG epoched\';
    base_output_path = 'C:\Users\sdabiri\OneDrive - Georgia Institute of Technology\Dr. Sederberg MaTRIX Lab\Research Paper\';
else
    error('Unknown path.');
end

% 1. Load Subject 14 P3 file
in_dir = fullfile(input_path, 'P3');
set_files = dir(fullfile(in_dir, '*.set'));
names_sorted = sort(cellstr({set_files.name})');
s14_file = names_sorted{14};

fprintf('Loading Subject 14: %s\n', s14_file);
EEG = pop_loadset('filename', s14_file, 'filepath', in_dir);

% 2. Run Extended Infomax ICA on Subject 14
fprintf('Running Infomax ICA on Subject 14 (All P3 Trials)...\n');
EEG = pop_runica(EEG, 'extended', 1, 'interupt', 'off');

% 3. Classify Components (Using EEGLAB ICLabel if installed, or inspect topoplots)
try
    EEG = iclabel(EEG);
    % Automatically flag ICs where Eye, Muscle, or Channel Noise probability > 0.70
    artifact_classes = EEG.etc.ic_classification.ICLabel.classifications;
    % Column 2: Muscle, Column 3: Eye, Column 5: Channel Noise
    bad_ics = find(artifact_classes(:, 2) > 0.70 | artifact_classes(:, 3) > 0.70 | artifact_classes(:, 5) > 0.70);
catch
    % Fallback: Pop up GUI to inspect component topographies
    pop_selectcomps(EEG, 1:min(12, size(EEG.icaweights, 1)));
    prompt = 'Enter bad IC indices to reject as a vector (e.g. [1 2]): ';
    bad_ics = input(prompt);
end

fprintf('Rejecting Artifact ICs: %s\n', mat2str(bad_ics));

% 4. Subtract Artifact Components
EEG_clean = pop_subcomp(EEG, bad_ics, 0);

% --- Ensure all EEGLAB headers (srate, times, pnts, xmin, xmax) are rebuilt ---
if ~isfield(EEG_clean, 'srate') || isempty(EEG_clean.srate) || EEG_clean.srate == 0
    EEG_clean.srate = EEG.srate;
end
if ~isfield(EEG_clean, 'xmin') || isempty(EEG_clean.xmin)
    EEG_clean.xmin = EEG.xmin;
    EEG_clean.xmax = EEG.xmax;
end

EEG_clean = eeg_checkset(EEG_clean); % Forces EEGLAB to rebuild EEG.times from srate and xmin

% 5. Save Cleaned Dataset to Disk
clean_save_dir = fullfile(input_path, 'P3_cleaned');
if ~exist(clean_save_dir, 'dir'), mkdir(clean_save_dir); end
pop_saveset(EEG_clean, 'filename', s14_file, 'filepath', clean_save_dir);

% 6. Quantitative Verification: Check Baseline Noise Ratio After Cleaning
time_s = linspace(-1.0, 2.5, EEG.pnts);
iBase  = time_s >= -1.0 & time_s < -0.8;

var_raw_before = mean(var(mean(EEG.data(:, iBase, :), 3), 0, 2));
var_raw_after  = mean(var(mean(EEG_clean.data(:, iBase, :), 3), 0, 2));

fprintf('\n--- VERIFICATION OF NOISE REDUCTION ---\n');
fprintf('Pre-cue Baseline Variance BEFORE ICA: %.4f uV^2\n', var_raw_before);
fprintf('Pre-cue Baseline Variance AFTER  ICA: %.4f uV^2\n', var_raw_after);
fprintf('Noise Reduction: %.1f%%\n', 100 * (var_raw_before - var_raw_after) / var_raw_before);