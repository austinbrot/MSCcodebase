function batch_MSC_analyses_BIDS_mod(current_MSCname_str)
% MODIFIED: Script to run MSC analyses, loading precomputed data.
% Expects a single string argument: the subject ID (e.g., 'MSC01')

if nargin < 1
    error('Usage: batch_MSC_analyses_BIDS_mod('current_MSCname_string')');
end

MSCname = strtrim(current_MSCname_str);
if isempty(MSCname)
    error('Received empty subject name string.');
end

fprintf('Starting batch MSC analysis for subject: %s\n', MSCname);

% All old logic for MSCnums_all, getenv('SUBJECT_INDEX_FROM_ARRAY') is removed.
% The script now processes ONLY the subject name passed as an argument.

edgethresh = .5;
xdist = 30;
thresholds = [.003 .004 .005:.005:.05];

%% Analyses to run
run_vertexwise_infomap = 1;
run_parcellation = 1;
make_parcel_corrmats = 1;
make_parcel_distmats = 1;
run_parcel_infomap = 1;
run_spring_embedding = 0;

home_dir = getenv('HOME');
oak_dir = getenv('OAK');
scratch_dir = getenv('SCRATCH');

if isempty(home_dir) || isempty(oak_dir) || isempty(scratch_dir)
    error('Required environment variables (HOME, OAK, SCRATCH) are not set.');
end

out_dir = fullfile(scratch_dir, '/MSCcodebase/test'); % Main output for this script's results

% NEW: Define path to precomputed data
precomputed_data_main_dir = fullfile(scratch_dir, 'MSCcodebase', 'precomputed_vertexwise_data');

MSC_dir = fullfile(oak_dir, '/inprocess/MSC/ds000224');
derivatives_dir = fullfile(oak_dir, '/inprocess/MSC/ds000224-derivatives-new');
surface_pipeine_dir = [derivatives_dir '/xcp_d']; % Corrected typo pipeline -> pipeline

surface_dist_dir = fullfile(oak_dir, '/inprocess/MSC/ds000224-derivatives/surface_pipeline');

sessions = {'01', '03', '05', '07', '09'};
% sessions = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10'};

% The main loop `for MSCnum = MSCnums` is no longer needed as we process one subject.
% The variable `MSCname` is now the definitive subject identifier from the argument.
    
fprintf('\nProcessing specific subject: %s\n', MSCname);

% NEW: Define paths for precomputed files for this subject
subject_precomp_dir = fullfile(precomputed_data_main_dir, MSCname);
alldata_file = fullfile(subject_precomp_dir, [MSCname '_alldata.mat']);
ciftitemplate_file = fullfile(subject_precomp_dir, [MSCname '_cifti_template.mat']);
vertex_corrmat_file = fullfile(subject_precomp_dir, [MSCname '_vertexwise_corrmat.mat']);

%% Define output folders for this subject's analysis results
infomap_outfolder = [out_dir '/infomap/' MSCname '_infomap_p003_p005_p05'];
parcellation_outfolder = [out_dir '/parcels'];
springembed_outfolder = [out_dir '/spring_embed'];
parcelinfomap_outfolder = [parcellation_outfolder '/' MSCname '_parcels_LR_infomap_p003_p05'];
surfdir = fullfile(derivatives_dir, 'fmriprep', ['sub-' MSCname], 'anat');

dmatname = fullfile(oak_dir, 'inprocess', 'MSC', 'fslr_distmat.mat');
    
% MODIFIED: Paths to ciftifiles still needed for parcel_creator_cifti template and potentially other steps
ciftifiles = cell(length(sessions),1);
rest_dir = fullfile(derivatives_dir, 'xcp_d', ['sub-' MSCname]);
for s = 1:length(sessions)
    ses_dir = fullfile(rest_dir, ['ses-func' sessions{s}], 'func');
    ciftifiles{s} = fullfile(ses_dir, ['sub-' MSCname '_ses-func' sessions{s} '_task-rest_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii']);
end
    
ciftidata_template_path = ciftifiles{1}; % Used as template for parcel_creator_cifti
parcellation_file = [parcellation_outfolder '/' MSCname '_parcels_LR.dtseries.nii'];

if ~exist(infomap_outfolder, 'dir'), mkdir(infomap_outfolder); end
if ~exist(springembed_outfolder, 'dir'), mkdir(springembed_outfolder); end
if ~exist(parcellation_outfolder, 'dir'), mkdir(parcellation_outfolder); end
if ~exist([parcellation_outfolder '/' MSCname '/'], 'dir'), mkdir([parcellation_outfolder '/' MSCname '/']); end
if ~exist(parcelinfomap_outfolder, 'dir'), mkdir(parcelinfomap_outfolder); end


% NEW: Load CIFTI metadata template (needed for writing outputs)
fprintf('Loading CIFTI metadata template from: %s\n', ciftitemplate_file);
if exist(ciftitemplate_file, 'file')
    load(ciftitemplate_file, 'cifti_metadata_template');
    data = cifti_metadata_template; % Use this as the base for cifti structures
else
    error('Precomputed CIFTI template file not found: %s', ciftitemplate_file);
end

%% vertexwise infomap
if run_vertexwise_infomap
    fprintf('Running vertex-wise Infomap for %s...\n', MSCname);
    cd(infomap_outfolder);
    
    % MODIFIED: Load precomputed vertex-wise corrmat
    fprintf('Loading precomputed vertex-wise corrmat from: %s\n', vertex_corrmat_file);
    if exist(vertex_corrmat_file, 'file')
        load(vertex_corrmat_file, 'corrmat'); % Loads 'corrmat' variable
        fprintf('Loaded corrmat with dimensions: %s\n', mat2str(size(corrmat)));
    else
        error('Precomputed vertex-wise corrmat file not found: %s', vertex_corrmat_file);
    end
    
    if ~isfield(data, 'brainstructure') || isempty(data.brainstructure)
         fprintf('Brainstructure not in template, attempting to load from %s\n', ciftifiles{1});
         if exist(ciftifiles{1}, 'file')
            temp_cifti_for_struct = ft_read_cifti_mod(ciftifiles{1});
            data.brainstructure = temp_cifti_for_struct.brainstructure;
            clear temp_cifti_for_struct;
         else
            error('Fallback CIFTI file for brainstructure not found: %s', ciftifiles{1});
         end
    end
    structure_indices = data.brainstructure;
    structure_indices = structure_indices(structure_indices>0);
    structure_indices = (structure_indices > 2) +1;
    
    Run_Infomap_2015(corrmat, dmatname, xdist, thresholds, 0, infomap_outfolder, 12, structure_indices);
    clear corrmat % Clear after use
    
    communities = modify_clrfile('simplify','rawassn.txt',400);
    regularized = regularize(communities);
    
    data_to_write = data;
    data_to_write.data = regularized;
    ft_write_cifti_mod([MSCname '_rawassn_minsize400_regularized'], data_to_write);
    try movefile('rawassn_minsize400_regularized.dtseries.nii',[MSCname '_rawassn_minsize400_regularized.dtseries.nii']); catch; end
    consensus_maker_knowncolors([MSCname '_rawassn_minsize400_regularized.dtseries.nii'],[],[],1);
    make_block_diagram([MSCname '_rawassn_minsize400_regularized_allcolumns_recolored.dscalar.nii'],thresholds);
    cifti_to_border_v2([MSCname '_rawassn_minsize400_regularized_recolored.dscalar.nii'],1,1,'default');
    fprintf('Finished vertex-wise Infomap for %s.\n', MSCname);
end

%% parcellation
if run_parcellation
    fprintf('Running parcellation for %s...\n', MSCname);
    cd([parcellation_outfolder '/' MSCname '/'])
    
    fprintf('Loading precomputed vertex-wise corrmat for parcellation from: %s\n', vertex_corrmat_file);
    if exist(vertex_corrmat_file, 'file')
        load(vertex_corrmat_file, 'corrmat');
    else
        error('Precomputed vertex-wise corrmat file not found: %s', vertex_corrmat_file);
    end
    
    data_for_parcellation = data;
    data_for_parcellation.data = corrmat;
    clear corrmat;
    
    fprintf('Parcellating surface for %s...\n', MSCname);
    surface_parcellation_singlesub_BIDS(MSCname, data_for_parcellation, surfdir, 100, 0, [parcellation_outfolder '/' MSCname]);
    clear data_for_parcellation;
    
    parcel_creator_cifti('corrofcorr_allgrad_LR_subcort_smooth2.55_wateredge_avg.dtseries.nii',[MSCname '_parcels'],edgethresh,ciftidata_template_path)
    movefile([parcellation_outfolder '/' MSCname '/' MSCname '_parcels_edgethresh_' num2str(edgethresh) '.dtseries.nii'],parcellation_file);
    try delete([parcellation_outfolder '/' MSCname '/corrofcorr_allgrad_LR_subcort_smooth2.55.dtseries.nii']); catch; end % MODIFIED: path to delete
    
    fprintf('Finished parcellation for %s!\n', MSCname);
end

%% parcel corrmats
if make_parcel_corrmats
    fprintf('Making parcel corrmats for %s...\n', MSCname);
    cd(parcellation_outfolder); % Change to general parcel outfolder

    parcels = ft_read_cifti_mod(parcellation_file);
    parcelIDs = unique(parcels.data); parcelIDs(parcelIDs<1) = [];
    
    fprintf('Loading precomputed alldata from: %s\n', alldata_file);
    if exist(alldata_file, 'file')
        load(alldata_file, 'alldata');
         fprintf('Loaded alldata with dimensions: %s\n', mat2str(size(alldata)));
    else
        error('Precomputed alldata file not found: %s', alldata_file);
    end
    
    tcs = zeros(size(alldata,2),length(parcelIDs));
    for IDnum = 1:length(parcelIDs)
        tcs(:,IDnum) = mean(alldata(parcels.data==parcelIDs(IDnum),:),1);
    end
    clear alldata;
    corrmat = paircorr_mod(tcs);
    corrmat(isnan(corrmat)) = 0;
    corrmat = FisherTransform(corrmat);
    
    save([parcellation_outfolder '/' MSCname '_parcel_corrmat.mat'],'corrmat');
    fprintf('Finished making parcel corrmats for %s.\n', MSCname);
end


%% parcel distmats
if make_parcel_distmats
    fprintf('Making parcel distmats for %s...\n', MSCname);
    cd(parcellation_outfolder); % Change to general parcel outfolder

    distances = smartload(dmatname);
    
    parcels = ft_read_cifti_mod(parcellation_file);
    parcelIDs = unique(parcels.data); parcelIDs(parcelIDs<1) = [];
    
    parcel_centroids = zeros(length(parcelIDs),1);
    for IDnum = 1:length(parcelIDs)
        parcelinds = find(parcels.data==parcelIDs(IDnum));
        within_parcel_distances = distances(parcelinds,parcelinds);
        [~,centroidind] = min(sum(within_parcel_distances,2));
        parcel_centroids(IDnum) = parcelinds(centroidind);
    end
    parcel_distances = distances(parcel_centroids,parcel_centroids);
    save([parcellation_outfolder '/' MSCname '_parcel_distances_xhemlarge.mat'],'parcel_distances');
    clear distances;
    fprintf('Finished making parcel distmats for %s.\n', MSCname);
end

%% parcel infomap
if run_parcel_infomap
    fprintf('Running parcel infomap for %s...\n', MSCname);
    cd(parcelinfomap_outfolder)

    load([parcellation_outfolder '/' MSCname '_parcel_corrmat.mat'], 'corrmat'); 
    load([parcellation_outfolder '/' MSCname '_parcel_distances_xhemlarge.mat'], 'distances');
    
    Run_Infomap_nopar(corrmat, distances, xdist, thresholds, 0, parcelinfomap_outfolder);
    communities = modify_clrfile('simplify','rawassn.txt',4);
    communities = load('rawassn_minsize4.txt');
    regularized = regularize(communities);
    
    parcels = ft_read_cifti_mod(parcellation_file);
    parcelIDs = unique(parcels.data); parcelIDs(parcelIDs<1) = [];
    
    regularized_out = data;
    regularized_out.data = zeros(size(parcels.data,1),size(regularized,2));
    for IDnum = 1:length(parcelIDs)
        regularized_out.data(parcels.data==parcelIDs(IDnum),:) = repmat(regularized(IDnum,:),nnz(parcels.data==parcelIDs(IDnum)),1);
    end
    ft_write_cifti_mod('rawassn_minsize4_regularized.dtseries.nii',regularized_out);
    consensus_maker_knowncolors('rawassn_minsize4_regularized.dtseries.nii',[],[infomap_outfolder '/' MSCname '_rawassn_minsize400_regularized_recolored.dscalar.nii'],1,[],[],parcellation_file);
    fprintf('Finished parcel infomap for %s.\n', MSCname);
end


%% spring embedding
if run_spring_embedding
    fprintf('Running spring embedding for %s...\n', MSCname);
    cd(springembed_outfolder)

    load([parcellation_outfolder '/' MSCname '_parcel_corrmat.mat'], 'corrmat');
    load([parcellation_outfolder '/' MSCname '_parcel_distances_xhemlarge.mat'], 'distances');
    consensus = load([parcelinfomap_outfolder '/rawassn_minsize4_regularized_recolored.txt']);
    spring_embedding_func_easy_crossthresh(corrmat,consensus,1,25,distances,xdist,thresholds,[MSCname '_spring_embed']);
    close all;
    fprintf('Finished spring embedding for %s.\n', MSCname);
end

fprintf('\nFinished all analyses for subject: %s\n', MSCname);
fprintf('Batch MSC analyses (modified) complete for %s.\n', MSCname);

end % End of function 