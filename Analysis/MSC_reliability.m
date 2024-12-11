
warning off
run_fromscratch = true;

MSCnums = [1];
totalMSCnums = 1;

thresholds = [.003 .004 .005:.005:.05];
sessions = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10'};

home_dir = getenv('HOME');
addpath(genpath(fullfile(home_dir, '2019_03_03_BCT')))
oak_dir = getenv('OAK');
scratch_dir = getenv('SCRATCH');
out_dir = fullfile(scratch_dir, 'MSCcodebase', 'results');
MSC_dir = fullfile(oak_dir, 'inprocess', 'MSC', 'ds000224');
derivatives_dir = fullfile(oak_dir, 'inprocess', 'MSC', 'ds000224-derivatives');
surface_pipeine_dir = fullfile(derivatives_dir, 'surface_pipeline');

convergence_out_dir = fullfile(out_dir, 'convergence');
mkdir(convergence_out_dir);
cd(convergence_out_dir);

% NOTE: The following data quantity parameters have been changed for faster
% testing. threshold and session have been modified above as well.
% The original values are commented out and should be replaced for
% the full analysis. The following files have also been modified:
% - run_infomap_on_pajekfile.m

xdistance = 30;
% iterations = 1000;
iterations = 100;
splithalf_quant = 70;
datalength_totest = [2.5 5 10 : 10 : 100];
TR = 2.4;
xdist = 30;


splithalf_quant_frames = round(splithalf_quant * 60 / TR);
datalength_totest_frames = round(datalength_totest * 60 / TR);
%%
if isempty(gcp('nocreate'))
    pool = parpool(12);
end

if run_fromscratch

    corrmat_similarity = zeros(iterations,length(datalength_totest),totalMSCnums);
    community_dice = zeros(iterations,length(datalength_totest),totalMSCnums);
    % PC_similarity = zeros(iterations,length(datalength_totest),totalMSCnums);
    GEff_delta = zeros(iterations,length(datalength_totest),totalMSCnums);
    Modularity_delta = zeros(iterations,length(datalength_totest),totalMSCnums);
    RichClub_delta = zeros(iterations,length(datalength_totest),totalMSCnums);
    GEff_all = zeros(iterations,length(datalength_totest),totalMSCnums);
    Modularity_all = zeros(iterations,length(datalength_totest),totalMSCnums);

    % corrmats_setaside = cell(iterations,length(datalength_totest),max(MSCnums));
    % corrmats_test = cell(iterations,length(datalength_totest),max(MSCnums));
    % communities_setaside = cell(iterations,length(datalength_totest),max(MSCnums));
    % communities_test = cell(iterations,length(datalength_totest),max(MSCnums));

else
    load([convergence_out_dir '/similarity_metrics.mat'])
    %load([convergence_out_dir '/corrmats_and_communities.mat'])
end
    

for MSCnum = MSCnums
    
    %if any(any(corrmat_similarity(:,:,MSCnum)==0))
    
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    
    ciftiparcelsdir = fullfile(out_dir, 'parcels');
    parcels_LR = [ciftiparcelsdir '/' MSCname '_parcels_LR.dtseries.nii'];
    parcels_struct = ft_read_cifti_mod(parcels_LR);
    parcels = parcels_struct.data;
    parcelIDs = unique(parcels); parcelIDs(parcelIDs<1) = [];
    parcel_distances = smartload([ciftiparcelsdir '/' MSCname '_parcel_distances_xhemlarge.mat']);
    alldata_communities = load(fullfile(ciftiparcelsdir, [ MSCname '_parcels_LR_infomap_p003_p05'], 'rawassn_minsize4_regularized_recolored.txt'));
    % PC_all{MSCnum} = zeros(length(parcelIDs),iterations,length(datalength_totest));
    
    cifti_all = [];
    tmask_all = [];
    sessnum = [];
    withinsess_ind = [];
    numframes = zeros(length(sessions),1);
    rest_dir = fullfile(surface_pipeine_dir, ['sub-' MSCname], 'processed_restingstate_timecourses');
    for s = 1 : length(sessions)
        ciftifiles{s} = fullfile(rest_dir, ['ses-func' sessions{s}], 'cifti', ['sub-' MSCname '_ses-func' sessions{s} '_task-rest_bold_32k_fsLR.dtseries.nii']);
        tmask_file = fullfile(rest_dir, ['ses-func' sessions{s}], 'cifti', ['sub-' MSCname '_ses-func' sessions{s} '_task-rest_bold_32k_fsLR_tmask.txt']);
        tmask = load(tmask_file);
        temp = ft_read_cifti_mod(ciftifiles{s});
        temp.data = temp.data(:,logical(tmask));
        tcs{s} = zeros(nnz(tmask),length(parcelIDs));
        for IDnum = 1:length(parcelIDs)
            tcs{s}(:,IDnum) = mean(temp.data(parcels==parcelIDs(IDnum),:),1);
        end
        numframes(s) = nnz(tmask);
        
    end
    
    
    numsessions = length(sessions);
    
    
    
    prevstring = [];
    
    for amountnum = 1:length(datalength_totest_frames)
        for iter = 1:iterations
            
            string = [MSCname ': iteration ' num2str(iter) ', ' num2str(datalength_totest(amountnum)) ' minutes of data'];
            disp(string)
            
            iteration_failed = true;
            
            while iteration_failed
                
                randorder = randperm(length(sessions));
                
                half_num_sessions = ceil(length(randorder)*1/2);
                
                %get set-aside data
                set_aside_sessions = randorder(1:half_num_sessions);
                set_aside_session_lengths = numframes(set_aside_sessions);
                [set_aside_session_lengths,sortorder] = sort(set_aside_session_lengths,'ascend');
                set_aside_sessions = set_aside_sessions(sortorder);
                
                set_aside_tcs = zeros(0,length(parcelIDs));
                
                for sesscount = 1:length(set_aside_sessions);
                    sessnum = set_aside_sessions(sesscount);
                    left_tograb = splithalf_quant_frames - size(set_aside_tcs,1);
                    thissess_length = size(tcs{sessnum},1);
                    amount_tograb_fromsess = min([thissess_length (ceil(left_tograb ./ (length(set_aside_sessions) - sesscount + 1)))]);
                    
                    wiggle_room_insess = thissess_length - amount_tograb_fromsess;
                    startpos = randi((wiggle_room_insess+1),1);
                    set_aside_tcs((end+1):(end+amount_tograb_fromsess),:) = tcs{sessnum}(startpos : (amount_tograb_fromsess + startpos - 1),:);
                end
                
                if size(set_aside_tcs,1) >= splithalf_quant_frames;
                    iteration_failed = false;
                end
            end
            
            
            %get test data
            test_sessions = setdiff([1:length(sessions)],set_aside_sessions);
            test_session_lengths = numframes(test_sessions);
            [test_session_lengths,sortorder] = sort(test_session_lengths,'ascend');
            test_sessions = test_sessions(sortorder);
            
            test_tcs = zeros(0,length(parcelIDs));
            
            for sesscount = 1:length(test_sessions);
                sessnum = test_sessions(sesscount);
                left_tograb = datalength_totest_frames(amountnum) - size(test_tcs,1);
                thissess_length = size(tcs{sessnum},1);
                amount_tograb_fromsess = min([thissess_length (ceil(left_tograb ./ (length(test_sessions) - sesscount + 1)))]);
                
                wiggle_room_insess = thissess_length - amount_tograb_fromsess;
                startpos = randi((wiggle_room_insess+1),1);
                test_tcs((end+1):(end+amount_tograb_fromsess),:) = tcs{sessnum}(startpos : (amount_tograb_fromsess + startpos - 1),:);
            end
            
            if size(test_tcs,1) < datalength_totest_frames(amountnum);
                iteration_failed = true;
            end
            
            
            if iteration_failed
                corrmat_similarity(iter,amountnum,MSCnum) = NaN;
                community_dice(iter,amountnum,MSCnum) = NaN;
                % PC_similarity(iter,amountnum,MSCnum) = NaN;
                GEff_delta(iter,amountnum,MSCnum) = NaN;
                Modularity_delta(iter,amountnum,MSCnum) = NaN;
                % PC_all{MSCnum}(:,iter,amountnum) = NaN;
                GEff_all(iter,amountnum,MSCnum) = NaN;
                Modularity_all(iter,amountnum,MSCnum) = NaN;
                
            else
                
                
                
                
                
                this_set_aside_outfolder = [convergence_out_dir '/setaside_infomap'];
                mkdir(this_set_aside_outfolder)
                this_test_outfolder = [convergence_out_dir '/test_infomap'];
                mkdir(this_test_outfolder)
                
                
                
                
                cd(this_set_aside_outfolder)
                
                set_aside_corrmat = paircorr_mod(set_aside_tcs);
                set_aside_corrmat(isnan(set_aside_corrmat)) = 0;
                set_aside_corrmat = FisherTransform(set_aside_corrmat);
                
                
                Run_Infomap_alreadyparpool(set_aside_corrmat, parcel_distances, xdistance, thresholds, 0, this_set_aside_outfolder)
                
                communities = modify_clrfile('simplify','rawassn.txt',4);
                regularized = regularize(communities);
                consensus_maker_knowncolors_textonly(regularized,alldata_communities,'rawassn_minsize4_regularized')
                set_aside_communities = load('rawassn_minsize4_regularized_recolored.txt');
                % [~, set_aside_meanPCs, ~] = PC_calc(set_aside_corrmat,parcel_distances < xdistance,communities,thresholds,ones(size(communities,1),1));
                
                globalEffs = zeros(length(thresholds),1);
                modularity = zeros(length(thresholds),1);
                richclub = zeros(length(thresholds),1);
                sorted_corr_vecs = sort(set_aside_corrmat(triu(true(size(set_aside_corrmat)),1)),'descend');
                for thresh = 1:length(thresholds)
                    rthresh = sorted_corr_vecs(round(length(sorted_corr_vecs) .* thresholds(thresh)));
                    corrmat_bin = set_aside_corrmat >= rthresh;
                    globalEffs(thresh) = efficiency_bin(corrmat_bin);
                    
                    [modularity(thresh), ~, ~, ~] = calc_modularity_TL(communities(:,thresh),corrmat_bin);
                    
                    
                end
                set_aside_GE = mean(globalEffs);
                set_aside_mod = mean(modularity);
                
                
                
                cd(this_test_outfolder)
                
                test_corrmat = paircorr_mod(test_tcs);
                test_corrmat(isnan(test_corrmat)) = 0;
                test_corrmat = FisherTransform(test_corrmat);
                
                Run_Infomap_alreadyparpool(test_corrmat, parcel_distances, xdistance, thresholds, 0, this_test_outfolder)
                communities = modify_clrfile('simplify','rawassn.txt',4);
                regularized = regularize(communities);
                consensus_maker_knowncolors_textonly(regularized,alldata_communities,'rawassn_minsize4_regularized')
                test_communities = load('rawassn_minsize4_regularized_recolored.txt');
                % [~, test_meanPCs, ~] = PC_calc(test_corrmat,parcel_distances < xdistance,communities,thresholds,ones(size(communities,1),1));
                
                % PC_all{MSCnum}(:,iter,amountnum) = test_meanPCs;
                
                globalEffs = zeros(length(thresholds),1);
                modularity = zeros(length(thresholds),1);
                richclub = zeros(length(thresholds),1);
                sorted_corr_vecs = sort(test_corrmat(triu(true(size(test_corrmat)),1)),'descend');
                for thresh = 1:length(thresholds)
                    rthresh = sorted_corr_vecs(round(length(sorted_corr_vecs) .* thresholds(thresh)));
                    corrmat_bin = test_corrmat >= rthresh;
                    globalEffs(thresh) = efficiency_bin(corrmat_bin);
                    
                    [modularity(thresh), ~, ~, ~] = calc_modularity_TL(communities(:,thresh),corrmat_bin);
                    
                    
                end
                test_GE = mean(globalEffs);
                test_mod = mean(modularity);
                
                GEff_all(iter,amountnum,MSCnum) = test_GE;
                Modularity_all(iter,amountnum,MSCnum) = test_mod;
                
                
                
                indmat = triu(true(size(test_corrmat)),1);
                
                corrmat_similarity(iter,amountnum,MSCnum) = paircorr_mod(set_aside_corrmat(indmat),test_corrmat(indmat));
                community_dice(iter,amountnum,MSCnum) = 2 * nnz(set_aside_communities==test_communities) ./ (numel(set_aside_communities) + numel(test_communities));
                % PC_naninds = isnan(set_aside_meanPCs) | isnan(test_meanPCs);
                % PC_similarity(iter,amountnum,MSCnum) = corr(set_aside_meanPCs(~PC_naninds),test_meanPCs(~PC_naninds));
                GEff_delta(iter,amountnum,MSCnum) = abs((test_GE - set_aside_GE)./set_aside_GE);
                Modularity_delta(iter,amountnum,MSCnum) = abs((test_mod - set_aside_mod)./set_aside_mod);
                
                
            end
            
            
            
            
        end
        
        save([convergence_out_dir '/similarity_metrics.mat'],'corrmat_similarity','community_dice','GEff_delta','Modularity_delta','RichClub_delta','GEff_all','Modularity_all')
    end
    
    disp(' ')
    
end

    
delete(gcp('nocreate'))

