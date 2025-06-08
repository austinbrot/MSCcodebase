function precompute_vertexwise_data(subject_list_csv_string)
% Script: precompute_vertexwise_data.m
% Purpose: Precompute alldata and vertex-wise correlation matrices for MSC subjects
% to reduce memory load in subsequent analysis steps.
% Accepts a comma-separated string of subject IDs (e.g., 'MSC01,MSC02,MSC03')

if nargin < 1
    error('Usage: precompute_vertexwise_data(''subject_id_list_csv'')');
end

fprintf('Starting precomputation of vertex-wise data...\n');
fprintf('Received subject list: %s\n', subject_list_csv_string);

%% Parse subject list
subject_ids_cell_array = strsplit(subject_list_csv_string, ',');
if isempty(subject_ids_cell_array) || all(cellfun('isempty', subject_ids_cell_array))
    error('Parsed subject ID list is empty. Check input CSV string.');
end
fprintf('Found %d subjects to process.\n', length(subject_ids_cell_array));

%% Define paths (environmental variables should be set by calling script)
home_dir = getenv('HOME');
oak_dir = getenv('OAK');
scratch_dir = getenv('SCRATCH');

if isempty(home_dir) || isempty(oak_dir) || isempty(scratch_dir)
    error('Required environment variables (HOME, OAK, SCRATCH) are not set.');
end

% Define output directory for precomputed files
precomputed_data_main_dir = fullfile(scratch_dir, 'MSCcodebase', 'precomputed_vertexwise_data');
if ~exist(precomputed_data_main_dir, 'dir')
    mkdir(precomputed_data_main_dir);
    fprintf('Created directory: %s\n', precomputed_data_main_dir);
end

MSC_dir = fullfile(oak_dir, 'data', 'MSC', 'ds000224'); % Base BIDS directory
derivatives_dir = fullfile(oak_dir, 'data', 'MSC', 'ds000224-derivatives'); % Derivatives directory

% Sessions to process for each subject
sessions = {'01', '03', '05', '07', '09'};
% sessions = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10'}; % Example for 10 sessions

%% Loop through subjects for precomputation
for i = 1:length(subject_ids_cell_array)
    MSCname = strtrim(subject_ids_cell_array{i}); % Use subject ID directly
    if isempty(MSCname)
        fprintf('Skipping empty subject ID at index %d.\n', i);
        continue;
    end
    fprintf('\nProcessing subject: %s\n', MSCname);

    subject_precomp_dir = fullfile(precomputed_data_main_dir, MSCname);
    if ~exist(subject_precomp_dir, 'dir')
        mkdir(subject_precomp_dir);
        fprintf('Created subject precomputation directory: %s\n', subject_precomp_dir);
    end

    alldata_file = fullfile(subject_precomp_dir, [MSCname '_alldata.mat']);
    ciftitemplate_file = fullfile(subject_precomp_dir, [MSCname '_cifti_template.mat']);
    vertex_corrmat_file = fullfile(subject_precomp_dir, [MSCname '_vertexwise_corrmat.mat']);

    %% Aggregate data from all sessions
    ciftifiles = cell(length(sessions),1);
    rest_dir = fullfile(derivatives_dir, 'xcpd-0.10.7', ['sub-' MSCname]);
    for s = 1:length(sessions)
        ses_dir = fullfile(rest_dir, ['ses-func' sessions{s}], 'func');
        ciftifiles{s} = fullfile(ses_dir, ['sub-' MSCname '_ses-func' sessions{s} '_task-rest_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii']);
    end

    fprintf('Loading and concatenating session data for %s...\n', MSCname);
    alldata = []; % Initialize alldata
    cifti_metadata_template = []; % Initialize template

    for s = 1:length(sessions)
        fprintf('Reading: %s\n', ciftifiles{s});
        if ~exist(ciftifiles{s}, 'file')
            warning('File not found: %s. Skipping this session for %s.', ciftifiles{s}, MSCname);
            continue;
        end
        data_struct = ft_read_cifti_mod(ciftifiles{s});
        session_data = data_struct.data;
        display(size(session_data));
        if isempty(alldata) % Check if alldata is empty for the first valid session
            alldata = session_data;
            cifti_metadata_template = data_struct; % Store the first one as template
            cifti_metadata_template.data = []; % Clear data from template
        else
            alldata = [alldata session_data];
        end
        clear session_data data_struct; % Clear to save memory
    end
    
    if isempty(alldata)
        warning('No data loaded for subject %s. Skipping precomputation for this subject.', MSCname);
        continue;
    end

    fprintf('Saving alldata for %s...\n', MSCname);
    save(alldata_file, 'alldata', '-v7.3');
    fprintf('Dimensions of alldata: %s\n', mat2str(size(alldata)));
    
    fprintf('Saving CIFTI metadata template for %s...\n', MSCname);
    save(ciftitemplate_file, 'cifti_metadata_template', '-v7.3');
    
    %% Compute and save vertex-wise correlation matrix
    fprintf('Computing vertex-wise corrmat for %s...\n', MSCname);
    corrmat = paircorr_mod(alldata'); % Transpose alldata to get (timepoints x vertices)
    clear alldata; % Clear alldata as soon as it's not needed
    
    corrmat(isnan(corrmat)) = 0;
    corrmat = FisherTransform(corrmat);
    
    fprintf('Saving vertex-wise corrmat for %s...\n', MSCname);
    save(vertex_corrmat_file, 'corrmat', '-v7.3');
    fprintf('Dimensions of corrmat: %s\n', mat2str(size(corrmat)));
    clear corrmat; % Clear corrmat after saving

    fprintf('Finished precomputation for %s.\n', MSCname);
end

fprintf('\nAll precomputations finished.\n');
end % End of function 