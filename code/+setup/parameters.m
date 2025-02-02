function [params, all_names] = parameters(runopts)
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    import solver.DHACalibrator
    
    %% SET SHARED PARAMETERS
    shared_params = struct();
    shared_params.r = 0.01;
    shared_params.Bequests = false;

    % Statistics from the 2019 SCF
    scf = setup.scf2019struct();
    shared_params.numeraire_in_dollars = scf.annual_earnings;
    shared_params.abars = [1000, 2000, 3000]  / scf.annual_earnings;

    mean_totw_target = scf.mean_totw / scf.annual_earnings;

    % Income processes
    quarterly_b_params = shared_params;
    quarterly_b_params.IncomeProcess = 'input/income_quarterly_b_contyT.mat';
    quarterly_b_params.sd_logyT = sqrt(0.6376);
    quarterly_b_params.lambdaT = 0.25;

    quarterly_c_params = shared_params;
    quarterly_c_params.IncomeProcess = 'input/income_quarterly_c_contyT.mat';
    quarterly_c_params.sd_logyT = sqrt(1.6243);
    quarterly_c_params.lambdaT = 0.0727;

    quarterly_a_params = shared_params;
    quarterly_a_params.sd_logyT = sqrt(0.2087);
    quarterly_a_params.sd_logyP = sqrt(0.01080);
    quarterly_a_params.rho_logyP = 0.9881;
    quarterly_a_params.lambdaT = 1;

    annual_params = shared_params;
    annual_params.sd_logyT = sqrt(0.0494);
    annual_params.sd_logyP = sqrt(0.0422);
    annual_params.rho_logyP = 0.9525;
    annual_params.lambdaT = 1;
    
    %% SET CALIBRATIONS
    % Ignored if not calibrating
    calibrations = struct();
    calibrations.variables = {'beta0'};
    % calibrations.target_names = {'median_a'};
    % calibrations.target_values = [scf.median_totw];
    calibrations.target_names = {'mean_a'};
    calibrations.target_values = [mean_totw_target];

    for ii = 2:999
        calibrations(ii) = calibrations(1);
    end

    %% PARAMETERIZATIONS

    %----------------------------------------------------------------------
    % BASELINES
    %----------------------------------------------------------------------
    
    % Annual
    params(1) = setup.Params(1, 'Annual', annual_params);
    params(1).beta0 = 0.984108034755346;
    params(1).group = [-1];
    params(1).colnums = [2];
    params(1).descr = 'Baseline';
    params(end).tex_header = 'Annual';
     
    % Quarterly
    params(end+1) = setup.Params(4, 'Quarterly', quarterly_b_params);
    params(end).beta0 = 0.984363510593659;
    params(end).group =     [1 2 3 4 5 6 7 8 -1 -2 -4];
    params(end).colnums =   [1 1 1 1 1 1 1 1  1  1  1];
    params(end).descr = 'Baseline';
    params(end).tex_header = 'Quarterly';

    %----------------------------------------------------------------------
    % PART 2, DIFFERENT ASSUMPTIONS
    %----------------------------------------------------------------------
    ifreq = 4;

    % Liquid wealth calibration, mean assets = 2.25
    name = sprintf('E[a] = %g', 2.25);
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).descr = 'Calibration to liquid wealth, E[a] = 2.25';
    params(end).tex_header = 'E[a]';
    params(end).tex_header_values.value = 2.25;
    n = numel(params);
    calibrations(n).target_names = {'mean_a'};
    calibrations(n).target_values = [2.25];
    
    % Total wealth calibration, mean assets = 9.4
    name = sprintf('E[a] = %g', 9.4);
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).group = [1];
    params(end).colnums = [2];
    params(end).descr = {'Calibration to total wealth, E[a] = 9.4'};
    params(end).tex_header = 'E[a]';
    params(end).tex_header_values.value = 9.4;
    params(end).betaH0 = 2e-3;
    n = numel(params);
    calibrations(n).target_names = {'mean_a'};
    calibrations(n).target_values = [9.4];
    
    % Liquid wealth calibration, median assets = 0.05, 0.5, 1.0
    for mw = [0.05, 0.5, 1.0]
        name = sprintf('median(a) = %g', mw);
        params(end+1) = setup.Params(4, name, quarterly_b_params);
        params(end).descr = {sprintf('Calibration to liquid wealth, median(a) = %g', mw)};
        params(end).tex_header = 'Median(a)';
        params(end).tex_header_values.value = mw;
        n = numel(params);
        calibrations(n).target_names = {'median_a'};
        calibrations(n).target_values = [mw];
    end

    % Other asset/HtM targets
    for mw = [0.046, 1.54]
        name = sprintf('median(a) = %g', mw);
        params(end+1) = setup.Params(4, name, quarterly_b_params);
        params(end).descr = {sprintf('Calibration to liquid wealth, median(a) = %g', mw)};
        params(end).tex_header = 'Median(a)';
        params(end).tex_header_values.value = mw;
        n = numel(params);
        calibrations(n).target_names = {'median_a'};
        calibrations(n).target_values = [mw];

        if mw == 1.54
            params(end).group = [1];
            params(end).colnums = [3];
        else
            params(end).group = [1];
            params(end).colnums = [5];
        end
    end

    name = sprintf('E[a] = %g', 0.5617);
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).group = [1];
    params(end).colnums = [4];
    params(end).descr = {'Calibration to total wealth, E[a] = 0.5617'};
    params(end).tex_header = 'E[a]';
    params(end).tex_header_values.value = 0.5617;
    params(end).betaH0 = 2e-3;
    n = numel(params);
    calibrations(n).target_names = {'mean_a'};
    calibrations(n).target_values = [0.5617];

    name = sprintf('HtM = %g', 0.142);
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).group = [1];
    params(end).colnums = [6];
    params(end).descr = {'Calibration to PHtM, HtM = 0.142'};
    params(end).tex_header = 'HtM';
    params(end).tex_header_values.value = 0.142;
    params(end).betaH0 = 2e-3;
    n = numel(params);
    calibrations(n).target_names = {'a_lt_ysixth'};
    calibrations(n).target_values = [0.142];

    name = sprintf('HtM = %g', 0.409);
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).descr = {'Calibration to PHtM, HtM = 0.409'};
    params(end).tex_header = 'HtM';
    params(end).tex_header_values.value = 0.409;
    params(end).betaH0 = 2e-3;
    n = numel(params);
    calibrations(n).target_names = {'a_lt_ysixth'};
    calibrations(n).target_values = [0.409];


    % no death
    name = 'No Death';
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).dieprob = 0;
    params(end).beta0 = 0.975363510593659;
    params(end).tex_header = 'No Death';
    params(end).group = [-2];
    params(end).colnums = [2];

    % % no bequests
    % name = 'No Bequests';
    % params(end+1) = setup.Params(4, name, quarterly_b_params);
    % params(end).group = {'Q1b'};
    % params(end).Bequests = 0;
    % params(end).tex_header = 'No Bequests';

    % With bequests
    name = 'With Bequests';
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).tex_header = 'With Bequests';
    params(end).group = [-2];
    params(end).colnums = [3];

    % perfect annuities
    name = 'Annuities';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).annuities = true;
    params(end).betaH0 = - 5e-3;
    params(end).tex_header = 'Annuities';
    params(end).group = [-2];
    params(end).colnums = [4];

    % different interest rates
    for ii = [0, 5]
        name = sprintf('r = %g%% p.a.', ii);
        params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
        params(end).r = ii/100;
        params(end).descr = name;

        params(end).tex_header = 'r';
        params(end).tex_header_values.r = sprintf('%g', ii/100);

        if ii == 0
            params(end).group = [2];
            params(end).colnums = [2];
        elseif ii == 5
            params(end).group = [2];
            params(end).colnums = [3];
        end
    end
    
    % interest rate heterogeneity
    name = 'Permanent r het, r in {-1,1,3} p.a.';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).r = [-1, 1, 3] / 100;
    params(end).betaH0 = -1e-4;
    params(end).beta0 = 0.973149481985717;
    params(end).descr = 'r in {-1, 1, 3}';
    params(end).tex_header = 'r';
    params(end).tex_header_values.r = '-0.01, 0.01, 0.03';
    params(end).group = [3];
    params(end).colnums = [6];
    
    name = 'Permanent r het, r in {-3,1,5} p.a.';
    params(end+1) = setup.Params(ifreq,name, quarterly_b_params);
    params(end).r = [-3, 1, 5] / 100;
    params(end).betaH0 = 1e-3;
    params(end).beta0 = 0.960885729527277;
    params(end).descr = 'r in {-3,1,5}';
    params(end).tex_header = 'r';
    params(end).tex_header_values.value = '-0.03, 0.01, 0.05';
    params(end).group = [3];
    params(end).colnums = [7];


%         % different tax rates
%         for itax = [0.05, 0.1, 0.15, 0.25]
%             name = [lfreq ' LabTax' num2str(itax)];
%             params(end+1) = setup.Params(ifreq,name,'');
%             params(end).labtaxlow = itax;
%         end


%         % bequest curvature
%         for bcurv = [0.1 0.5 1 2 5]
%             name = [lfreq ' BeqWt0.02 BeqLux0.01 BeqCurv' num2str(bcurv)];
%             params(end+1) = setup.Params(ifreq,name,IncomeProcess);
%             params(end).bequest_weight = 0.02;
%             params(end).bequest_luxury = 0.01;
%             params(end).bequest_curv   = bcurv;target_value
%         end
   
    % fixed beta heterogeneity
    for ibw = [0.001, 0.005, 0.01]
        name = sprintf('Beta5, pSwitch0, pSpacing%g', ibw);
        params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
        params(end).nbeta = 5;
        params(end).betawidth = ibw;
        params(end).prob_zswitch = 0;
        params(end).beta0 = 0.956194383870642;
        params(end).descr = sprintf('p = %g, spacing = %g', 0, ibw);
        params(end).tex_header = '$\beta$ het.';
        params(end).tex_header_values.pswitch = 0;
        params(end).tex_header_values.spacing = ibw;
        
        if ibw == 0.005
            params(end).betaH0 = 7e-4;
        elseif ibw == 0.01
            params(end).betaH0 = 1e-3;
        end

        if ibw == 0.005
            params(end).group = [3];
            params(end).colnums = [2];
        elseif ibw == 0.01
            params(end).group = [3];
            params(end).colnums = [3];
        end
    end

    % fixed beta het w/normal
    ibw = 0.01;
    name = sprintf('Beta5, pSwitch0 (normal), pSpacing%g', ibw);
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).nbeta = 5;
    params(end).betawidth = ibw;
    params(end).prob_zswitch = 0;
    params(end).beta0 = 0.956194383870642;
    params(end).zdist_forced = [0.1, 0.2, 0.4, 0.2, 0.1];
    params(end).descr = sprintf('p = %g, normal, spacing = %g', 0, ibw);
    params(end).tex_header = '$\beta$ het.';
    params(end).tex_header_values.pswitch = 0;
    params(end).tex_header_values.spacing = ibw;
    params(end).betaH0 = 1e-2;

    % random beta heterogeneity
    for ibw = [0.01]
        for bs = [1/50, 1/10]
            name = sprintf('Beta5, pSwitch%g, pSpacing%g', bs, ibw);
            params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
            params(end).nbeta = 5;
            params(end).betawidth = ibw;
            params(end).prob_zswitch = bs;
            params(end).descr = sprintf('p = %g, spacing = %g', bs, ibw);
            params(end).tex_header = '$\beta$ het.';
            params(end).tex_header_values.pswitch = bs;
            params(end).tex_header_values.spacing = ibw;

            if bs == 1 / 50
                params(end).betaH0 = -1.2e-2;
            else
                params(end).betaH0 = -1e-3;
            end

            if (ibw == 0.01) && (bs == 1/50)
                params(end).group = [3];
                params(end).colnums = [4];
            elseif (ibw == 0.01) && (bs == 1/10)
                params(end).group = [3];
                params(end).colnums = [5];
            end
        end
    end

    % Different risk aversion coeffs
    for ira = [0.5, 2, 6]
        name = sprintf('CRRA = %g', ira);
        params(end+1) = setup.Params(4, name, quarterly_b_params);
        params(end).risk_aver = ira;
        if (ifreq==4 && ira==4) || ira==6
            params(end).betaL = 0.5;
        end

        if ira == 6
            params(end).betaH0 = -1e-3;
        end
        params(end).descr = 'CRRA';
        params(end).tex_header = 'CRRA';
        params(end).tex_header_values.value = ira;

        if ira == 0.5
            params(end).group = [2];
            params(end).colnums = [4];
        elseif ira == 6
            params(end).group = [2];
            params(end).colnums = [5];
        end
    end
    
    % CRRA with IES heterogeneity
    name = 'CRRA w/IES betw exp(-1), exp(1)';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).risk_aver = 1./ exp([-1 -0.5 0 0.5 1]);
    if params(end).freq == 4
        params(end).betaH0 =  - 2e-3;
    end
    params(end).descr = 'RA = exp(1), ..., exp(-1), IES = exp(-1), ..., exp(1)';
    params(end).tex_header = 'CRRA';
    params(end).tex_header_values.riskaver = 'exp(-1), ..., exp(1)';
    
    name = 'CRRA w/IES betw exp(-2), exp(2)';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).risk_aver = 1 ./ exp([-2 -1 0 1 2]);
    if params(end).freq == 4
        params(end).betaH0 = -1e-3;
    end
    params(end).descr = 'RA = exp(2), ..., exp(-2), IES = exp(-2), ..., exp(2)';
    params(end).tex_header = 'CRRA';
    params(end).tex_header_values.riskaver = 'exp(-2), ..., exp(2)';

    name = 'CRRA w/IES betw exp(-3), exp(3)';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).risk_aver = 1 ./ exp([-3 -1.5 0 1.5 3]);
    if params(end).freq == 4
        params(end).betaH0 = -1e-3;
    end
    params(end).descr = 'RA = exp(3), ..., exp(-3), IES = exp(-3), ..., exp(3)';
    params(end).tex_header = 'CRRA';
    params(end).tex_header_values.riskaver = 'exp(-3), ..., exp(3)';

    % epstein-zin, quarterly
    ras = [0.5 8  1    1]; % 8
    ies = [1   1  0.25 2]; % 2
    for iez = 1:4
        ra_i = ras(iez);
        ies_i = ies(iez);
        name = sprintf('EZ, ra%g, ies%g', ra_i, ies_i);
        params(end+1) = setup.Params(4, name, quarterly_b_params);
        params(end).risk_aver = ra_i;
        params(end).invies = 1 / ies_i;
        params(end).EpsteinZin = true;
        params(end).descr = sprintf('RA = %g, IES = %g', ra_i, ies_i);
        params(end).tex_header = 'EZ';
        params(end).tex_header_values.riskaver = ra_i;
        params(end).tex_header_values.ies = ies_i;

        switch iez
            case 1
                params(end).betaH0 = - 2e-3;
            case 2
                params(end).betaH0 = - 3e-3;
            case 3
                params(end).betaH0 = - 3e-3;
            case 4
                params(end).betaH0 = - 1e-3;
            case 5
                params(end).betaH0 = 3e-2;
                params(end).beta0 = 0.999868 ^ 4;
        end

        params(end).group = [-4];
        params(end).colnums = [iez + 1];
    end

    % EZ with IES heterogeneity
    name = 'EZ w/ IES betw exp(-1), exp(1)';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).invies = 1 ./ exp([-1 -0.5 0 0.5 1]);
    params(end).EpsteinZin = true;
    if (ifreq == 4)
        params(end).betaH0 = - 3e-3;
    end
    params(end).descr = 'RA = 1, IES = exp(-1), ..., exp(1)';
    params(end).tex_header = 'EZ';
    params(end).tex_header_values.riskaver = 1;
    params(end).tex_header_values.ies = 'exp(-1), ..., exp(1)';
    params(end).group = [4];
    params(end).colnums = [2];
    
    name = 'EZ w/ IES betw exp(-2), exp(2)';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).invies = 1 ./ exp([-2 -1 0 1 2]);
    params(end).EpsteinZin = true;
    if (ifreq == 4)
        params(end).betaH0 = - 3e-3;
    end
    params(end).descr = 'RA = 1, IES = exp(-2), ..., exp(2)';
    params(end).tex_header = 'EZ';
    params(end).tex_header_values.riskaver = 1;
    params(end).tex_header_values.ies = 'exp(-2), ..., exp(2)';
    params(end).group = [4];
    params(end).colnums = [3];

    name = 'EZ w/ IES betw exp(-3), exp(3)';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).invies = 1 ./ exp([-3 -1.5 0 1.5 3]);
    params(end).EpsteinZin = true;
    if (ifreq == 4)
        params(end).betaH0 = - 3e-3;
    end
    params(end).descr = 'RA = 1, IES = exp(-3), ..., exp(3)';
    params(end).tex_header = 'EZ';
    params(end).tex_header_values.riskaver = 1;
    params(end).tex_header_values.ies = 'exp(-3), ..., exp(3)';
    params(end).group = [4];
    params(end).colnums = [4];

    % EZ with risk aversion heterogeneity
    name = 'EZ w/ RA betw exp(-2), exp(2)';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).invies = 1;
    params(end).risk_aver = exp([-2 -1 0 1 2]);
    params(end).EpsteinZin = 1;
    params(end).betaH0 = - 1e-3;
    params(end).betaL = 0.96;
    params(end).beta0 = 0.999349 ^ 4;
    params(end).descr = 'RA = exp(-2), ..., exp(2), IES = 1';
    params(end).tex_header = 'EZ';
    params(end).tex_header_values.ies = 1;
    params(end).tex_header_values.riskaver = 'exp(-2), ..., exp(2)';

    name = 'EZ w/ RA betw exp(-3), exp(3)';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).invies = 1;
    params(end).risk_aver = exp([-3 -1.5 0 1.5 3]);
    params(end).EpsteinZin = 1;
    params(end).betaH0 = - 1e-3;
    params(end).betaL = 0.96;
    params(end).beta0 = 0.999349 ^ 4;
    params(end).descr = 'RA = exp(-3), ..., exp(3), IES = 1';
    params(end).tex_header = 'EZ';
    params(end).tex_header_values.ies = 1;
    params(end).tex_header_values.riskaver = 'exp(-3), ..., exp(3)';

    % temptation
    for tempt = [0.01 0.05 0.07]
        name = sprintf('Temptation = %g', tempt);
        params(end+1) = setup.Params(4, name, quarterly_b_params);
        params(end).temptation = tempt;
        params(end).tex_header = 'Temptation';
        params(end).tex_header_values.tempt = sprintf('%g', tempt);

        if tempt == 0.01
            params(end).group = [4];
            params(end).colnums = [5];
        elseif tempt == 0.05
            params(end).group = [4];
            params(end).colnums = [6];
            params(end).betaH0 = -2e-5;
        elseif tempt == 0.07
            params(end).beta0 = 1;
            params(end).betaH0 = 2e-3;
        end

    end

    name = 'Temptation, uniform in {0, 0.05, 0.1}';
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).temptation = [0 0.05 0.1];
    params(end).beta0 = 0.999083 ^ 4;
    params(end).betaH0 = -4e-5;
    params(end).descr = 'Temptation in {0, 0.05, 0.1}';
    params(end).tex_header = 'Temptation';
    params(end).tex_header_values.tempt = '0, 0.05, 0.1';
    params(end).group = [4];
    params(end).colnums = [7];


    %----------------------------------------------------------------------
    % PART 3a, ANNUAL MODEL
    %----------------------------------------------------------------------
    
    % i
    name = 'Annual, no yT';
    params(end+1) = setup.Params(1, name, annual_params);
    params(end).beta0 = 0.99;
    params(end).nyT = 1;
    params(end).sd_logyT = 0;
    params(end).lambdaT = 0;
    params(end).descr = 'No trans shocks';
    params(end).tex_header = 'Annual (i)';
    params(end).tex_header_values.descr = 'No trans shocks';

    % iv
    name = 'Annual, Carrol';
    params(end+1) = setup.Params(1, name, annual_params);
    params(end).rho_logyP = 0.999;
    params(end).sd_logyP = sqrt(0.015);
    params(end).sd_logyT = sqrt(0.01);
    params(end).descr = 'Carrol process';
    params(end).tex_header = 'Annual (iv)';
    params(end).tex_header_values.descr = 'Carrol process';
    
%     % v
%     params(end+1) = setup.Params(1,'A a(v) HighPersNotReEst','');
%     params(end).rho_logyP = 0.99;
%     
%     % vi
%     params(end+1) = setup.Params(1,'A a(vi) LowPersNotReEst','');
%     params(end).rho_logyP = 0.9;
% 
%     % vii
%     params(end+1) = setup.Params(1,'A a(vii) HighPersReEst','');
%     params(end).rho_logyP = 0.99;
%     params(end).sd_logyP = sqrt(0.0088);
%     params(end).sd_logyT = sqrt(0.0667);
    
    % viii
    name = 'Annual, high persistence';
    params(end+1) = setup.Params(1, name, shared_params);
    params(end).rho_logyP = 0.995;
    params(end).sd_logyP = sqrt(0.0043);
    params(end).sd_logyT = sqrt(0.0688);
    params(end).lambdaT = 1;
    params(end).descr = 'High persistence';
    params(end).tex_header = 'Annual (viii)';
    params(end).tex_header_values.descr = 'High persistence';
%     
%     % ix
%     params(end+1) = Params(1,'A a(ix) HighPersNoTransReEst','');
%     params(end).rho_logyP = 0.99;
%     params(end).sd_logyP = sqrt(0.0088);
%     params(end).nyT = 1;
%     params(end).sd_logyT = sqrt(0);
    
    % x
    name = 'Annual, high nyF = 5';
    params(end+1) = setup.Params(1, name, shared_params);
    params(end).rho_logyP = 0.9158;
    params(end).sd_logyP = sqrt(0.0445);
    params(end).sd_logyT = sqrt(0.0479);
    params(end).sd_logyF = sqrt(0.1801);
    params(end).lambdaT = 1;
    params(end).nyF = 5;
    params(end).descr = 'FE heterogeneity';
    params(end).tex_header = 'Annual (x)';
    params(end).tex_header_values.descr = 'FE heterogeneity';
    
    %----------------------------------------------------------------------
    % PART 3b, QUARTERLY MODEL
    %----------------------------------------------------------------------
    
    % i quarterly_a
    name = 'quarterly_a';
    params(end+1) = setup.Params(4, name, quarterly_a_params);
    params(end).beta0 = 0.984363510593659;
    params(end).descr = 'quart_a';
    params(end).tex_header = 'Quart (i)';
    params(end).tex_header_values.descr = 'quart$_a$';
    params(end).group = [-3];
    params(end).colnums = [1];
    
    % ii
    name = 'KMP';
    params(end+1) = setup.Params(4, name, shared_params);
    params(end).rho_logyP = 0.9879;
    params(end).sd_logyP = sqrt(0.0109);
    params(end).sd_logyT = sqrt(0.0494);
    params(end).descr = 'KMP';
    params(end).tex_header = 'Quart (ii)';
    params(end).tex_header_values.descr = 'KMP';
    params(end).group = [-3];
    params(end).colnums = [2];

    % iii quarterly_c
    name = 'quarterly_c';
    params(end+1) = setup.Params(4, name, quarterly_c_params);
    params(end).beta0 = 0.984363510593659;
    params(end).descr = 'quart_c';
    params(end).tex_header = 'Quart (iii)';
    params(end).tex_header_values.descr = 'quart$_c$';
    params(end).group = [-3];
    params(end).colnums = [3];
    
    % no transitory shocks
    name = 'no trans shocks';
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).beta0 = 0.984363510593659;
    params(end).nyT = 1;
    params(end).descr = 'no_trans_shocks';
    params(end).tex_header = 'Quart No Trans';
    params(end).tex_header_values.descr = 'quart no trans';
    params(end).group = [-3];
    params(end).colnums = [4];
    
    %----------------------------------------------------------------------
    % EXTRAS
    %----------------------------------------------------------------------
    ibw = 0.01;
    name = sprintf('Beta5, pSwitch0, 0.05-0.1-0.7-0.1-0.05, pSpacing%g', ibw);
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).nbeta = 5;
    params(end).betawidth = ibw;
    params(end).prob_zswitch = 0;
    params(end).beta0 = 0.956194383870642;
    params(end).zdist_forced = [0.05, 0.1, 0.7, 0.1, 0.05];
    params(end).descr = sprintf('p = %g, normal, spacing = %g', 0, ibw);
    params(end).tex_header = '$\beta$ het.';
    params(end).tex_header_values.pswitch = 0;
    params(end).tex_header_values.spacing = ibw;
    params(end).betaH0 = 1e-2;
    
    
    for betaLfixed = 0.2:0.2:0.8
        name = sprintf('betaL = %g, betaH calibrated', betaLfixed);
        params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
        params(end).nbeta = 2;
        params(end).prob_zswitch = 0;
        params(end).beta0 = betaLfixed;
        params(end).betaL = params(end).beta0;
        params(end).betaHvariable = betaLfixed + 0.002;
        params(end).zdist_forced = [0.15, 0.85];
        params(end).descr = name;
        params(end).tex_header = '$\beta$ het.';
        params(end).tex_header_values.pswitch = 0;
        params(end).tex_header_values.spacing = ibw;
        params(end).betaH0 = 1e-2;
        n = numel(params);
        calibrations(n).variables = {'betaHvariable'};
    end

    name = '3pt CRRA = [exp(-2), exp(0), exp(2)], p = [0.15, 0.7, 0.15]';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).risk_aver = exp([-2 0 2]);
    params(end).zdist_forced = [0.15, 0.7, 0.15];
    params(end).betaH0 = -1e-3;
    params(end).descr = name;
    params(end).tex_header = '$\beta$ het.';
    params(end).tex_header_values.pswitch = 0;    


    name = '2pt CRRA = [exp(-3), exp(0)], p = [0.15, 0.85]';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).risk_aver = [exp(-3), exp(0)];
    params(end).zdist_forced = [0.15, 0.85];
    params(end).betaH0 = -1e-3;
    params(end).descr = name;
    params(end).tex_header = '$\beta$ het.';
    params(end).tex_header_values.pswitch = 0;   

    name = 'Permanent r het, r in {-1,1,3} p.a., p = [0.15, 0.7, 0.15]';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).r = [-1, 1, 3] / 100;
    params(end).betaH0 = 1e-3;
    params(end).beta0 = 0.997149481985717;
    params(end).zdist_forced = [0.15, 0.7, 0.15];
    params(end).descr = 'r in {-1, 1, 3}';
    params(end).tex_header = 'r';
    params(end).tex_header_values.r = '-0.01, 0.01, 0.03';

    name = 'Permanent r het, r in {-1,3} p.a., p = [0.15, 0.85]';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).r = [-1, 3] / 100;
    params(end).betaH0 = 1e-3;
    params(end).beta0 = 0.973149481985717;
    params(end).zdist_forced = [0.15, 0.85];
    params(end).descr = 'r in {-1, 3}';
    params(end).tex_header = 'r';
    params(end).tex_header_values.r = '-0.01, 0.03';
    
    % Quarterly baseline quad_b = 1
    params(end+1) = setup.Params(4, 'Quarterly quad_b = 1', quarterly_b_params);
%     params(end).beta0 = 0.984363510593659;
    params(end).beta0 = 1.01;
    params(end).betaH0 = 2e-3;
%     params(end).group =     [1 2 3 4 5 6 7 8 -1 -2 -4];
%     params(end).colnums =   [1 1 1 1 1 1 1 1  1  1  1];
    params(end).descr = 'Baseline, b = 1';
    params(end).tex_header = 'Quarterly';
    params(end).quad_b = 1.0;
    
    % Quarterly baseline quad_b = 0.5
    params(end+1) = setup.Params(4, 'Quarterly quad_b = 0.5', quarterly_b_params);
%     params(end).beta0 = 0.984363510593659;
    params(end).beta0 = 1.01;
    params(end).betaH0 = 2e-3;
%     params(end).group =     [1 2 3 4 5 6 7 8 -1 -2 -4];
%     params(end).colnums =   [1 1 1 1 1 1 1 1  1  1  1];
    params(end).descr = 'Baseline, b = 0.5';
    params(end).tex_header = 'Quarterly';
    params(end).quad_b = 0.5;
    
    % Quarterly baseline quad_b = 0.25
    params(end+1) = setup.Params(4, 'Quarterly quad_b = 0.25', quarterly_b_params);
%     params(end).beta0 = 0.984363510593659;
    params(end).beta0 = 1.01;
    params(end).betaH0 = 2e-3;
%     params(end).group =     [1 2 3 4 5 6 7 8 -1 -2 -4];
%     params(end).colnums =   [1 1 1 1 1 1 1 1  1  1  1];
    params(end).descr = 'Baseline, b = 0.25';
    params(end).tex_header = 'Quarterly';
    params(end).quad_b = 0.25;
    
    % Quarterly baseline quad_b = 0.1
    params(end+1) = setup.Params(4, 'Quarterly quad_b = 0.1', quarterly_b_params);
%     params(end).beta0 = 0.984363510593659;
    params(end).beta0 = 1.01;
    params(end).betaH0 = 2e-3;
%     params(end).group =     [1 2 3 4 5 6 7 8 -1 -2 -4];
%     params(end).colnums =   [1 1 1 1 1 1 1 1  1  1  1];
    params(end).descr = 'Baseline, b = 0.1';
    params(end).tex_header = 'Quarterly';
    params(end).quad_b = 0.1;
    
    % Quarterly baseline quad_b = 0.05
    params(end+1) = setup.Params(4, 'Quarterly quad_b = 0.05', quarterly_b_params);
%     params(end).beta0 = 0.984363510593659;
    params(end).beta0 = 1.01;
    params(end).betaH0 = 2e-3;
%     params(end).group =     [1 2 3 4 5 6 7 8 -1 -2 -4];
%     params(end).colnums =   [1 1 1 1 1 1 1 1  1  1  1];
    params(end).descr = 'Baseline, b = 0.05';
    params(end).tex_header = 'Quarterly';
    params(end).quad_b = 0.05;
    
    % Quarterly baseline exp_a = 2
    params(end+1) = setup.Params(4, 'Quarterly exp_a = 2', quarterly_b_params);
%     params(end).beta0 = 0.984363510593659;
    params(end).beta0 = 1.01;
    params(end).betaH0 = 2e-3;
%     params(end).group =     [1 2 3 4 5 6 7 8 -1 -2 -4];
%     params(end).colnums =   [1 1 1 1 1 1 1 1  1  1  1];
    params(end).descr = 'Baseline, exp_a = 2';
    params(end).tex_header = 'Quarterly';
    params(end).exp_a = 2;
    
    % Quarterly baseline exp_a = 1
    params(end+1) = setup.Params(4, 'Quarterly exp_a = 1', quarterly_b_params);
%     params(end).beta0 = 0.984363510593659;
    params(end).beta0 = 1.01;
    params(end).betaH0 = 2e-3;
%     params(end).group =     [1 2 3 4 5 6 7 8 -1 -2 -4];
%     params(end).colnums =   [1 1 1 1 1 1 1 1  1  1  1];
    params(end).descr = 'Baseline, exp_a = 1';
    params(end).tex_header = 'Quarterly';
    params(end).exp_a = 1;
    
    % Quarterly baseline exp_a = 0.5
    params(end+1) = setup.Params(4, 'Quarterly exp_a = 0.5', quarterly_b_params);
%     params(end).beta0 = 0.984363510593659;
    params(end).beta0 = 1.01;
    params(end).betaH0 = 2e-3;
%     params(end).group =     [1 2 3 4 5 6 7 8 -1 -2 -4];
%     params(end).colnums =   [1 1 1 1 1 1 1 1  1  1  1];
    params(end).descr = 'Baseline, exp_a = 0.5';
    params(end).tex_header = 'Quarterly';
    params(end).exp_a = 0.5;
    
    % Quarterly baseline exp_a = 0.2
    params(end+1) = setup.Params(4, 'Quarterly exp_a = 0.2', quarterly_b_params);
%     params(end).beta0 = 0.984363510593659;
    params(end).beta0 = 1.01;
    params(end).betaH0 = 2e-3;
%     params(end).group =     [1 2 3 4 5 6 7 8 -1 -2 -4];
%     params(end).colnums =   [1 1 1 1 1 1 1 1  1  1  1];
    params(end).descr = 'Baseline, exp_a = 0.2';
    params(end).tex_header = 'Quarterly';
    params(end).exp_a = 0.2;
    
    % Quarterly baseline exp_a = 0.1
    params(end+1) = setup.Params(4, 'Quarterly exp_a = 0.1', quarterly_b_params);
%     params(end).beta0 = 0.984363510593659;
    params(end).beta0 = 1.01;
    params(end).betaH0 = 2e-3;
%     params(end).group =     [1 2 3 4 5 6 7 8 -1 -2 -4];
%     params(end).colnums =   [1 1 1 1 1 1 1 1  1  1  1];
    params(end).descr = 'Baseline, exp_a = 0.1';
    params(end).tex_header = 'Quarterly';
    params(end).exp_a = 0.1;
    
    % E[a] = 0.56, quad_b = 1
    name = sprintf('E[a] = %g, quad_b = 1', 0.5617);
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).group = [1];
    params(end).colnums = [4];
    params(end).descr = {'Calibration to total wealth, E[a] = 0.5617, b = 1'};
    params(end).tex_header = 'E[a]';
    params(end).tex_header_values.value = 0.5617;
    params(end).betaH0 = 2e-3;
    params(end).quad_b = 1.0;
    n = numel(params);
    calibrations(n).target_names = {'mean_a'};
    calibrations(n).target_values = [0.5617];
    
    % fixed beta het w/normal, quad_b = 1
    ibw = 0.01;
    name = sprintf('Beta5, pSwitch0 (normal), pSpacing%g', ibw);
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).nbeta = 5;
    params(end).betawidth = ibw;
    params(end).prob_zswitch = 0;
    params(end).beta0 = 0.956194383870642;
    params(end).quad_b = 1.0;
    params(end).zdist_forced = [0.1, 0.2, 0.4, 0.2, 0.1];
    params(end).descr = sprintf('p = %g, normal, spacing = %g, b = 1', 0, ibw);
    params(end).tex_header = '$\beta$ het.';
    params(end).tex_header_values.pswitch = 0;
    params(end).tex_header_values.spacing = ibw;
    params(end).betaH0 = 1e-2;
    
    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS, DO NOT CHANGE
    %----------------------------------------------------------------------
    params.set_index();

    % get list of all names
    all_names = cell2table({params.name}');
    
    % select by number if there is one, otherwise select by names,
    % otherwise use all
    if numel(runopts.number) == 1
        params = params(runopts.number);
    elseif numel(runopts.number) > 1
        error('runopts.number must have 1 or zero elements')
    else
        params = setup.Params.select_by_name(params, runopts.name_to_run);
    end

    params.set_run_parameters(runopts);
    params.make_adjustments();

    %----------------------------------------------------------------------
    % ATTACH CALIBRATOR
    %----------------------------------------------------------------------
    if params.calibrating
        calibration = calibrations(params.index);
        calibrator = DHACalibrator(params, calibration.variables,...
            calibration.target_names, calibration.target_values);
        heterogeneity = setup.Prefs_R_Heterogeneity(params);

        if isempty(params.betaHvariable) && (params.nbeta > 1) ...
                && isequal(heterogeneity.ztrans, eye(params.nbeta))
            % Adjust upper bound on beta if using fixed discount factor heterogeneity
            new_betaH = params.betaH - max(heterogeneity.betagrid0);
            params.set("betaH", new_betaH, true);
        end

        beta_bounds = [params.betaL, params.betaH];
        calibrator.set_param_bounds(beta_bounds);
        calibrator.set_handle(params);
        params.set("calibrator", calibrator, true);
    end
end
