classdef TexTables
    
    properties (Constant)
        table_includes;
    end

    methods (Static)
        function [n, labels, vars, decimals] = get_header_variables(tableno)
            decimals = [];

            if tableno == 1
                labels = {'Value'};
                vars = {'value'};
                decimals = [3];
            elseif tableno == 2
                labels = {'Value'};
                vars = {'value'};
                decimals = [2];
            elseif tableno == 3
                labels = {'Beta Spacing', 'r', 'Switching Probability'};
                vars = {'spacing', 'r', 'pswitch'};
            elseif tableno == 4
                labels = {'Risk Aversion', 'IES', 'Temptation'};
                vars = {'riskaver', 'ies', 'tempt'};
                decimals = [2, 3, 2];
            elseif tableno == -1
                labels = {};
                vars = {};
            elseif tableno == -2
                labels = {};
                vars = {};
            elseif tableno == -3
                labels = {'Description'};
                vars = {'descr'};
            elseif tableno == -4
                labels = {'Risk Aversion', 'IES'};
                vars = {'riskaver', 'ies'};
                decimals = [2, 2];
            end

            n = numel(labels);
            if isempty(decimals)
                decimals = ones(1, n);
            end
        end

        function save_baselines_tables(params_in, results, dirpath, varargin)
            for tableno = [1, 2]
	            for panelname = {'header', 'A', 'B', 'C', 'D'}
	            	if (tableno == 1)
	            		panelobj = tables.TexTables.table1panel(params_in, results, panelname{:}, varargin{:});
	            	elseif (tableno == 2)
	            		panelobj = tables.TexTables.table2panel(params_in, results, panelname{:}, varargin{:});
	            	end

	            	if strcmp(panelname{:}, 'header')
	            		panelfname = sprintf('table%d_header.xlsx', tableno);
	            	else
	            		panelfname = sprintf('table%d_panel%s.xlsx', tableno, panelname{:});
	            	end

	            	panelfpath = fullfile(dirpath, panelfname);
	            	writetable(panelobj, panelfpath, 'WriteRowNames', true);
	            end
	        end
        end

        function save_experiment_table(params_in, results, comparison_decomps, dirpath, tableno)
            header = tables.TexTables.experiment_table_header(params_in, results, tableno);
            tabcapnum = tablecnum(tableno);
            headerpath = fullfile(dirpath, sprintf('table%d_header.xlsx', tabcapnum));
            writetable(header, headerpath, 'WriteRowNames', true);

            for panelname = {'A', 'B', 'C', 'D'}
            	if ismember(panelname{:}, {'A', 'A2'})
            		panelobj = tables.TexTables.experiment_table_panel(...
            			params_in, comparison_decomps, panelname{:}, tableno);
            	else
            		panelobj = tables.TexTables.experiment_table_panel(...
            			params_in, results, panelname{:}, tableno);
            	end
            	panelfname = sprintf('table%d_panel%s.xlsx', tabcapnum, panelname{:});
            	panelfpath = fullfile(dirpath, panelfname);
            	writetable(panelobj, panelfpath, 'WriteRowNames', true);
            end
        end

        function table_out = table1panel(params_in, results, panel, varargin)
            switch panel
	            case 'header'
	            	get_stats = @(x) {
	            		x.stats.mpcs(5).quarterly
	                    x.stats.mpcs(5).annual
                        x.stats.mpcs(5).quarterly_htm_a_lt_1000
                        x.stats.mpcs(5).quarterly_mean_a
	                    x.stats.beta_A_effective
	                  };
	            case 'A'
		            get_stats = @(x) {...
		            	x.stats.mean_gross_y_annual
                        x.stats.std_log_gross_y_annual
                        x.stats.std_log_net_y_annual
	                };
	            case 'B'
	            	get_stats = @(x) {...
	            		x.stats.mean_a
                        x.stats.median_a
                        x.stats.sav0
                        x.stats.constrained{1}
                        x.stats.constrained_dollars{1}
                        x.stats.constrained_dollars{2}
                        x.stats.constrained_dollars{3}
                        x.stats.constrained_dollars{4}
                        x.stats.a_lt_ysixth
                        x.stats.a_lt_ytwelfth
                        x.stats.wpercentiles{1}
                        x.stats.wpercentiles{2}
                        x.stats.wpercentiles{4}
                        x.stats.wpercentiles{5}
                        x.stats.wpercentiles{7}
                        x.stats.wpercentiles{8}
                        x.stats.w_top10share
                        x.stats.w_top1share
                        x.stats.wgini
                    };
                case 'C'
                	get_stats = @(x) {...
                		x.stats.mpcs(4).annual
                        x.stats.mpcs(6).annual
                        x.stats.mpcs(4).quarterly
                        x.stats.mpcs(6).quarterly
                    };
                case 'D'
                	get_stats = @(x) {...
                		x.stats.mpcs(1).annual
                        x.stats.mpcs(2).annual
                        x.stats.mpcs(3).annual
                        x.stats.mpcs(1).quarterly
                        x.stats.mpcs(2).quarterly
                        x.stats.mpcs(3).quarterly
                    };
                otherwise
                	error("Invalid panel entry")
            end

            [stats_array, names_array] = stack_results(1, get_stats, params_in, results, varargin{:});
            table_out = make_table(stats_array, names_array);
        end

        function table_out = table2panel(params_in, results, panel, varargin)
            decomp_norisk_get_stats_fn = @(x, k)  {
            	x.stats.decomp_norisk.term2(k)
                x.stats.decomp_norisk.term3(k)
                x.stats.decomp_norisk.term4(k)
            };

            switch panel
	            case 'header'
	            	get_stats = @(x) {
	            		x.stats.mpcs(5).oneperiod
	            		x.stats.decomp_norisk.term1_pct
	            	};
	            case 'A'
	                get_stats = @(x) decomp_norisk_get_stats_fn(x, 1);
	            case 'B'
	            	get_stats = @(x) decomp_norisk_get_stats_fn(x, 2);
	           	case 'C'
	           		get_stats = @(x) decomp_norisk_get_stats_fn(x, 3);
	            case 'D'
	            	get_stats = @(x) {
	            		x.stats.decomp_RA.Em1_less_mRA
	                    x.stats.decomp_RA.term1
	                    x.stats.decomp_RA.term2
	                    x.stats.decomp_RA.term3
	                };
	           	otherwise
	           		error("Invalid panel selection")
           	end

            [stats_array, names_array] = stack_results(2, get_stats, params_in, results, varargin{:});
            table_out = make_table(stats_array, names_array);
        end

        function table_out = experiment_table_header(params_in, results, tableno)
            indices = filter_param_group(params_in, tableno);

            import statistics.Statistics.sfill

            [nhvars, hlabels, hvars, hdecimals] = tables.TexTables.get_header_variables(tableno);

            all_annual = true;
            for ip = indices
                if params_in(ip).freq == 4
                    all_annual = false;
                    break
                end
            end

            for ii = 1:numel(indices)
                ip = indices(ii);
                tex_vals = params_in(ip).tex_header_values;

                variable_values = {};
                for ih = 1:nhvars
                    variable_values{end+1} = sfill(tex_vals.(hvars{ih}), hlabels{ih}, hdecimals(ih));
                end

                if all_annual
                    statistics{ii} = {  results(ip).stats.mpcs(5).annual
                                        results(ip).stats.beta_A_effective
                                      };
                else
                    statistics{ii} = {  results(ip).stats.mpcs(5).quarterly
                                        results(ip).stats.mpcs(5).annual
                                        results(ip).stats.mpcs(5).quarterly_htm_a_lt_1000
                                        results(ip).stats.mpcs(5).quarterly_mean_a
                                        results(ip).stats.beta_A_effective
                                      };
                end
                statistics{ii} = [variable_values(:); statistics{ii}];
                names{ii} = params_in(ip).tex_header;
            end

            table_out = make_table(statistics, names, 'experiment', true);
        end

        function table_out = experiment_table_panel(params_in, variables, panel, tableno)
            indices = filter_param_group(params_in, tableno);
            all_annual = true;
            for ip = indices
                if params_in(ip).freq == 4
                    all_annual = false;
                    break
                end
            end

            switch panel
	            case 'A'
	            	get_stats = @(x) {
	            		x.Em1_less_Em0
                        x.term1
                        x.term2
                        x.term2a(1)
                        x.term2b(1)
                        x.term3
	            	};
	           	case 'A2'
	           		get_stats = @(x) {
	           			x.term1_pct
                        x.term2_pct
                        x.term2a_pct(1)
                        x.term2b_pct(1)
                        x.term3_pct
	           		};
	           	case 'B'
	           		get_stats = @(x) {
	           			x.stats.mean_a
                        x.stats.median_a
	                    x.stats.a_lt_ysixth
	                    x.stats.constrained_dollars{1}
	                    x.stats.constrained_dollars{2}
	                    x.stats.constrained_dollars{3}
	                    x.stats.constrained_dollars{4}
	                    x.stats.constrained_dollars{5}
	                    x.stats.w_top10share
	                    x.stats.w_top1share
	                    x.stats.wgini
	           		};
	           	case 'C'
                    if all_annual
                        get_stats = @(x) {
                            x.stats.mpcs(4).annual
                            x.stats.mpcs(6).annual
                        };
                    else
    	           		get_stats = @(x) {
    	           			x.stats.mpcs(4).quarterly
                            x.stats.mpcs(6).quarterly
                        };
                    end
                case 'D'
                    if all_annual
                        get_stats = @(x) {
                            x.stats.mpcs(1).annual
                            x.stats.mpcs(2).annual
                            x.stats.mpcs(3).annual
                        };
                    else
    	           		get_stats = @(x) {
    	           			x.stats.mpcs(1).quarterly
                            x.stats.mpcs(2).quarterly
                            x.stats.mpcs(3).quarterly
                        };
                    end
            end

            [stats_array, names_array] = stack_results(tableno, get_stats, params_in, variables, 'experiment', true);
            table_out = make_table(stats_array, names_array, 'experiment', true);
        end
    end 
end

function tnum = tablecnum(tableno)
    if tableno < 0
        tnum = sprintf('A%d', abs(tableno));
    else
        tnum = sprintf('%d', tableno);
    end
end

function indices = filter_param_group(params_in, tableno)
    ips = [];
    colnums = [];
    jj = 1;
    for ii = 1:numel(params_in)
        if ismember(tableno, params_in(ii).group)
            idx = find(tableno == params_in(ii).group);
            colnums = [colnums, params_in(ii).colnums(idx)];
            ips = [ips, ii];
        end
    end

    sorted_mat = sortrows([colnums(:), ips(:)]);
    indices = reshape(sorted_mat(:,2), 1, []);
end

function [stats_array, names_array] = stack_results(tableno, fn_handle, params_in, main_results, varargin)
    parser = inputParser;
    addOptional(parser, 'experiment', false);
    addOptional(parser, 'ctimeresults', []);
    parse(parser, varargin{:});
    experiment = parser.Results.experiment;
    ctimeresults = parser.Results.ctimeresults;

    indices = filter_param_group(params_in, tableno);
    n = numel(indices);
    stats_array = cell(n, 1);
    names_array = cell(n, 1);
    for ii = 1:n
        ip = indices(ii);

        idx = find(tableno == params_in(ip).group);
        colnum = params_in(ip).colnums(idx);
        stats_array{colnum} = fn_handle(main_results(ip));

        if experiment & ~isempty(params_in(ip).tex_header)
            names_array{colnum} = params_in(ip).tex_header;
        else
            names_array{colnum} = params_in(ip).name;
        end
    end

    if ~isempty(ctimeresults)
        if ismember(tableno, ctimeresults.p.group)
            idx = find(tableno == ctimeresults.p.group);
            colnum = ctimeresults.p.colnums(idx);
            stats_array{colnum} = fn_handle(ctimeresults);
            names_array{colnum} = 'Continuous Time';
        end
    end
end

function table_out = make_table(statistics, names, varargin)
    parser = inputParser;
    addOptional(parser, 'experiment', false);
    parse(parser, varargin{:});
    experiment = parser.Results.experiment;

    n = numel(names);
    for ii = 1:n
        vars{ii} = get_values(statistics{ii});

        if experiment
            varnames{ii} = char(names{ii} + string(sprintf('__v%d__', ii)));
        else
            varnames{ii} = names{ii};
        end
    end

    vars{n+1} = get_precision(statistics{1});
    varnames{n+1} = 'decimals';
    rownames = get_names(statistics{1});
    table_out = table(vars{:}, 'RowNames', rownames(:), 'VariableNames', varnames(:));
end

function values = get_values(entries)
    values = {};
    for ii = 1:numel(entries)
        values{ii} = entries{ii}.value;
    end
    values = values(:);
end

function decimals = get_precision(entries)
    decimals = {};
    for ii = 1:numel(entries)
        decimals{ii} = entries{ii}.decimals;
    end
    decimals = decimals(:);
end

function names = get_names(entries)
    names = {};
    for ii = 1:numel(entries)
        names{ii} = entries{ii}.tex_label;
    end
    names = names(:);
end