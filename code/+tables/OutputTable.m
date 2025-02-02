
function table_out = OutputTable(p, stats, comparison_decomps)
	if ~isempty(p.tex_header)
        name = p.tex_header;
    else
        name = p.name;
    end

	%% Intro panel
	descr = struct('value', p.descr, 'label', 'Description');
	intro_stats = {...
		descr
		stats.mpcs(5).quarterly
		stats.mpcs(5).quarterly_htm_biw
		stats.mpcs(5).annual
		stats.mpc_apc_corr{5}
		stats.beta_A
		stats.beta_A_effective
		stats.mean_gross_y_annual
		stats.std_log_gross_y_annual
		stats.std_log_net_y_annual
	};

	values = get_values(intro_stats);
	rownames = get_names(intro_stats);
	table_out = table(values, 'RowNames', rownames(:), 'VariableNames', {name});

	%% Wealth stats
	paneltitle = struct('value', NaN, 'label', '____Wealth Statistics');
	wealth_stats = {...
		paneltitle
		stats.mean_a
		stats.median_a
		stats.sav0
		stats.constrained_pct{1}
		stats.constrained_pct{2}
		stats.constrained_pct{3}
		stats.constrained_pct{4}
		stats.constrained_pct{5}
		stats.constrained_pct{6}
		stats.constrained_pct{7}
		stats.a_lt_ysixth
		stats.a_lt_ytwelfth
		stats.w_top10share
		stats.w_top1share
		stats.wgini
		stats.wpercentiles{1}
		stats.wpercentiles{2}
		stats.wpercentiles{3}
		stats.wpercentiles{5}
		stats.wpercentiles{7}
		stats.wpercentiles{8}
		stats.constrained_dollars{1}
		stats.constrained_dollars{2}
		stats.constrained_dollars{3}
		stats.constrained_dollars{4}
		stats.constrained_dollars{5}
	};

	values = get_values(wealth_stats);
	rownames = get_names(wealth_stats);
	new_rows = table(values, 'RowNames', rownames(:), 'VariableNames', {name});
	table_out = [table_out; new_rows];

	%% MPC Size Effects
	paneltitle = struct('value', NaN, 'label', '____MPC Size Effects');
	mpcsize_stats = {
		paneltitle
		stats.mpcs(4).quarterly
		stats.mpcs(6).quarterly
		stats.mpcs(4).annual
		stats.mpcs(6).annual
	};

	values = get_values(mpcsize_stats);
	rownames = get_names(mpcsize_stats);
	new_rows = table(values, 'RowNames', rownames(:), 'VariableNames', {name});
	table_out = [table_out; new_rows];

	%% MPC Sign Effects
	paneltitle = struct('value', NaN, 'label', '____MPC Sign Effects');
	mpcsign_stats = {
		paneltitle
		stats.mpcs(1).quarterly
		stats.mpcs(2).quarterly
		stats.mpcs(3).quarterly
		stats.mpcs(1).annual
		stats.mpcs(2).annual
		stats.mpcs(3).annual
	};

	values = get_values(mpcsign_stats);
	rownames = get_names(mpcsign_stats);
	new_rows = table(values, 'RowNames', rownames(:), 'VariableNames', {name});
	table_out = [table_out; new_rows];

	%% Intertemporal MPCs
	paneltitle = struct('value', NaN, 'label', '____MPC, Intertemporal');
	it_mpc_stats = {
		paneltitle
		stats.mpcs(5).avg_s_t{1,1}
		stats.mpcs(5).avg_s_t{1,2}
		stats.mpcs(5).avg_s_t{1,3}
		stats.mpcs(5).avg_s_t{1,4}
		stats.mpcs(5).avg_s_t{1,5}
	};

	values = get_values(it_mpc_stats);
	rownames = get_names(it_mpc_stats);
	new_rows = table(values, 'RowNames', rownames(:), 'VariableNames', {name});
	table_out = [table_out; new_rows];

	%% News MPCs
	paneltitle = struct('value', NaN, 'label', '____MPC, Out of News');
	mpcnews_stats = {
		paneltitle
		stats.mpcs(5).avg_s_t{2,1}
		stats.mpcs(5).avg_s_t{3,1}
		stats.mpcs(5).avg_s_t{4,1}
		stats.mpcs(5).avg_s_t{5,1}
	};

	values = get_values(mpcnews_stats);
	rownames = get_names(mpcnews_stats);
	new_rows = table(values, 'RowNames', rownames(:), 'VariableNames', {name});
	table_out = [table_out; new_rows];

	%% Decomp wrt RA model
	paneltitle = struct('value', NaN, 'label', '____Decomposition of E[MPC] - MPC_RA');
	decomp_ra_stats = {
		paneltitle
		stats.decomp_RA.Em1_less_mRA
		stats.decomp_RA.term1
        stats.decomp_RA.term2
        stats.decomp_RA.term3
    };

    values = get_values(decomp_ra_stats);
	rownames = get_names(decomp_ra_stats);
	new_rows = table(values, 'RowNames', rownames(:), 'VariableNames', {name});
	table_out = [table_out; new_rows];

	%% Decomp wrt no-risk model
	% Header
	paneltitle = struct('value', NaN, 'label',...
		'____Decomps of E[MPC] wrt RA and no inc risk, $500 shock');

	if p.freq == 1
		tmp = struct('value', stats.mpcs(5).annual.value, 'label', 'MPC (%)');
	else
		tmp = struct('value', stats.mpcs(5).quarterly.value, 'label', 'MPC (%)');
	end
	decomp_norisk_stats = {
		paneltitle
		tmp
		stats.decomp_norisk.term1_pct
    };

    values = get_values(decomp_norisk_stats);
	rownames = get_names(decomp_norisk_stats);
	new_rows = table(values, 'RowNames', rownames(:), 'VariableNames', {name});
	table_out = [table_out; new_rows];

	% HtM thresholds
	for ithresh = 1:numel(p.abars)
		threshold = p.abars(ithresh);
		panel_name = sprintf('__ _For HtM threshold %g', threshold);
		paneltitle = struct('value', NaN, 'label', panel_name);

		new_entries = {
			paneltitle
			stats.decomp_norisk.term2(ithresh)
			stats.decomp_norisk.term3(ithresh)
			stats.decomp_norisk.term4(ithresh)
		};

		values = get_values(new_entries);
		rownames = get_names(new_entries);
		new_rows = table(values, 'RowNames', rownames(:), 'VariableNames', {name});
		table_out = [table_out; new_rows];
	end

	%% Comparison decomps
	if nargin >= 3
		paneltitle = struct('value', NaN, 'label',...
		'____Decomps of E[MPC] wrt E[MPC_0], $500 shock');

		comparison_decomp_header = {
			paneltitle
			comparison_decomps.Em1_less_Em0
            comparison_decomps.term1
            comparison_decomps.term2
            comparison_decomps.term3
	    };

	    values = get_values(comparison_decomp_header);
		rownames = get_names(comparison_decomp_header);
		new_rows = table(values, 'RowNames', rownames(:), 'VariableNames', {name});
		table_out = [table_out; new_rows];

		for ithresh = 1:numel(p.abars)
			threshold = p.abars(ithresh);
			panel_name = sprintf('____Term 2 decomp for HtM threshold %g', threshold);
			paneltitle = struct('value', NaN, 'label', panel_name);

			new_entries = {
				paneltitle
                comparison_decomps.term2a(ithresh)
                comparison_decomps.term2b(ithresh)
			};

			values = get_values(new_entries);
			rownames = get_names(new_entries);
			new_rows = table(values, 'RowNames', rownames(:), 'VariableNames', {name});
			table_out = [table_out; new_rows];
		end
	end
end

function values = get_values(entries)
    values = {};
    for ii = 1:numel(entries)
        values{ii} = entries{ii}.value;
    end
    values = values(:);
end

function names = get_names(entries)
    names = {};
    for ii = 1:numel(entries)
        names{ii} = entries{ii}.label;
    end
    names = names(:);
end