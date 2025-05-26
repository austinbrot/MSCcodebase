file_new = '/scratch/users/abrotman/MSCcodebase/results-new/convergence/similarity_metrics.mat';
file_old = '/scratch/users/abrotman/MSCcodebase/results/convergence/similarity_metrics.mat';

% Load both datasets
data_new = load(file_new);
data_old = load(file_old);

% Get field names
names = fieldnames(data_new);
disp(names);

for j = 1:length(names)
    if ~strcmp(names{j}, 'PC_all')
        % Process new data
        eval([names{j} ' = reshape(data_new.' names{j} ', [size(data_new.' names{j} ',1), size(data_new.' names{j} ',2), 1]);']);
        eval([names{j} '(' names{j} ' == 0) = NaN;']);
        
        % Process old data
        eval([names{j} '_old = reshape(data_old.' names{j} ', [size(data_old.' names{j} ',1), size(data_old.' names{j} ',2), 1]);']);
        eval([names{j} '_old(' names{j} '_old == 0) = NaN;']);
    end
end

MSCnums = 1:1;
datalength_totest = [2.5 5 10:10:100];
outfolder = pwd;
colors_new = [0 0 1];  % Blue for "fmriprep + xcp-d"
colors_old = [1 0 0];  % Red for "MSC preprocessed"

%% Plotting
h = figure('Color', 'white', 'position', [1982 478 1352 804], 'DefaultAxesFontSize', 40);
hold on

legendnames = {'MSC preprocessed', 'fmriprep + xcp-d'};

% Loop for each metric
metrics = {'corrmat_similarity', 'community_dice', 'GEff_delta', 'Modularity_delta', 'GEff_all', 'Modularity_all'};
ylabels = {'Correlation to other half', 'Dice to other half', '% Difference from other half', '% Difference from other half', ...
           'Global Efficiency value', 'Modularity Value'};
titles = {'Connectivity Matrix', 'Network Assignment', 'Global Efficiency', 'Modularity', 'Global Efficiency', 'Modularity'};
ylimits = {[.1 .95], [.1 .95], [0 20], [0 20], [0 Inf], [0 Inf]};

for metric_idx = 1:length(metrics)
    metric = metrics{metric_idx};
    ylabel_text = ylabels{metric_idx};
    title_text = titles{metric_idx};
    ylim_values = ylimits{metric_idx};
    
    h = figure('Color', 'white', 'position', [1982 478 1352 804], 'DefaultAxesFontSize', 40);
    hold on
    
    % New Data
    for MSCnum = MSCnums
        eval(['meancorrmat_new = nanmean(' metric '(:,:,MSCnum), 1);']);
        eval(['meancorrmat_new(sum(isnan(' metric '(:,:,MSCnum)),1) > 900) = NaN;']);
        eval(['stdcorrmat_new = nanstd(' metric ', 0, 1);']);
        plot(datalength_totest, meancorrmat_new, 'Color', colors_new, 'LineWidth', 5);
    end

    % Old Data
    for MSCnum = MSCnums
        eval(['meancorrmat_old = nanmean(' metric '_old(:,:,MSCnum), 1);']);
        eval(['meancorrmat_old(sum(isnan(' metric '_old(:,:,MSCnum)),1) > 900) = NaN;']);
        eval(['stdcorrmat_old = nanstd(' metric '_old, 0, 1);']);
        plot(datalength_totest, meancorrmat_old, 'Color', colors_old, 'LineWidth', 5, 'LineStyle', '--');
    end
    
    % Set properties
    set(gca, 'FontSize', 40, 'FontWeight', 'bold', 'LineWidth', 3);
    title(gca, title_text);
    xlabel('Time (minutes)', 'FontWeight', 'bold', 'FontSize', 50);
    ylabel(ylabel_text, 'FontWeight', 'bold', 'FontSize', 50);
    ylim(ylim_values);
    legend(legendnames, 'Location', 'SouthEast');
    
    % Save the plot
    filename = [outfolder '/' title_text '.pdf'];
    try
        export_fig(gca, filename);
    catch
        savefig(gcf, [outfolder '/' title_text '.fig']);
    end
end
