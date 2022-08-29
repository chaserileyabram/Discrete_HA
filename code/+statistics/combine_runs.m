
% This script combines .mat files named variablesX.mat into
% an excel spreadsheet

clear

[~, currdir] = fileparts(pwd());
if ~strcmp(currdir, 'Discrete_HA')
    msg = 'The user must cd into the Discrete_HA directory';
    bad_dir = MException('Discrete_HA:master', msg);
    throw(bad_dir);
end

% Check if code is running locally or on the server
taskid = str2num(getenv('SLURM_ARRAY_TASK_ID'));
running_on_server = ~isempty(taskid);

if ~running_on_server
    taskid = 1;
end

outdir = fullfile('output', sprintf('tables%d', taskid));
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

addpath('code');

% Read .mat files into a cell array
for irun = 1:999
    fname = sprintf('variables%d.mat', irun);
    fpath = fullfile('output', fname);

    if ~exist(fpath, 'file')
        params{irun} = [];
        heterogeneity{irun} = [];
        stats{irun} = [];
        continue
    end

    S = load(fpath);
    params{irun} = S.Sparams;
    heterogeneity{irun} = S.heterogeneity;
    stats{irun} = S.results.stats;

    % Output tables
    xlxname = sprintf('table%d.xlsx', irun);
    savexlxpath = fullfile(outdir, xlxname);

    % Comparison decomps
    if params{irun}.freq == 1
        baseind = 1;
    else
        baseind = 2;
    end

    cdecomp = statistics.ComparisonDecomp(params{baseind}, params{irun},...
        stats{baseind}, stats{irun});

    mpcs0 = reshape(stats{baseind}.mpcs(5).mpcs_1_t{1},...
        params{baseind}.nx_DST, []);
    mpcs1 = reshape(stats{irun}.mpcs(5).mpcs_1_t{1},...
        params{irun}.nx_DST, []);
    cdecomp.perform_decompositions(mpcs0, mpcs1);

    table_out = tables.OutputTable(S.Sparams, S.results.stats, cdecomp.results);
    writetable(table_out, savexlxpath, 'WriteRowNames', true);
end

% Parameters table
params_cleaned = cell(1,1);
ip = 0;
for ii = 1:numel(params)
    if ~isempty(params{ii})
        ip = ip + 1;
        params_cleaned{ip} = params{ii};
    end
end
ptable = tables.param_tables(params_cleaned, heterogeneity);
panelfpath = fullfile(outdir, 'params_table.xlsx');
writetable(ptable, panelfpath, 'WriteRowNames', true);
