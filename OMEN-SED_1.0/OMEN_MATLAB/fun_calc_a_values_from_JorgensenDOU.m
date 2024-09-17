

% calls this one here:
% [res] = calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, Zinf, toc_load)


Db_Middelburg = false;  % we use our new fit, see /Documents/Projects/18_RECCAP2/OMEN/Bioturbation/Solan_ea_2019/make_empirical_fct.m
Zinf = 800;     % enough fo r

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% TODO: loop to calc a-values for uncertainty propagation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UNCERTAINTY VALUES

%% Restreppo SAR
SAR = 'Restreppo';
string_out = 'DOU_calc_Restrep_SA_TOCPOR';

toc_load = 'SA'
por_in = 'SA'
TOC_SA = struct2array(load('/home/domhu/Documents/GitHub/CoastalCarbonBurial/OMEN/CalcAvalues/output/Random_numbers_TOC_x5_240917.mat'));
POR_SA = struct2array(load('/home/domhu/Documents/GitHub/CoastalCarbonBurial/OMEN/CalcAvalues/output/Random_numbers_POR_x5_240917.mat'));
SA_exps = length(TOC_SA);
for k=1:SA_exps
    res_Restreppo_best = benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load, por_in, k, TOC_SA{k}, POR_SA{k});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE a-VALUES FOR BEST TOC & POR -- and also extreme BCs

% orginal toc from MArkus without back calculation
if false % used before, e.g., 240810
    
%% Restreppo SAR
SAR = 'Restreppo';
string_out = 'DOU_calc_Restrep';

toc_load = 'best'
por_in = 'por_best'
res_Restreppo_best = benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load, por_in);


SAR = 'Restreppo_low';
string_out = 'DOU_calc_Restrep_low';

toc_load = 'best'
por_in = 'por_best'
res_Restreppo_low_best = benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load, por_in);

SAR = 'Restreppo_high';
string_out = 'DOU_calc_Restrep_high';

toc_load = 'best'
por_in = 'por_best'
res_Restreppo_high_best = benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load, por_in);


%% Restreppo SAR
SAR = 'Restreppo';
string_out = 'DOU_calc_Restrep';

toc_load = 'best'
por_in = 'por_best'
res_Restreppo_best = benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load, por_in);

toc_load = 'best'
por_in = 'por_low'
res_Restreppo_best = benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load, por_in);

toc_load = 'best'
por_in = 'por_high'
res_Restreppo_best = benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load, por_in);


toc_load = 'low'
por_in = 'por_low'
res_Restreppo_low = benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load, por_in);

toc_load = 'high'
por_in = 'por_high'
res_Restreppo_high = benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load, por_in);



%% Burwicz SAR  -- only do best estimates here
SAR='Burwicz';
string_out = 'DOU_calc_Burwicz';

toc_load = 'best'
por_in = 'por_best'
res_Burwicz_best = benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load, por_in);


%% Middelburg SAR  -- only do best estimates here
SAR='Middelburg';
string_out = 'DOU_calc_Middelburg';

toc_load = 'best'
por_in = 'por_best'
res_Middelburg_best = benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load, por_in);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% OLD STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if false
    %% Burwicz SAR  -- only do best estimates here
    SAR='Burwicz';
    string_out = 'DOU_calc_Burwicz';
    
    % toc_load = 'best'
    % res_Burwicz_best = benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load);
    
    toc_load = 'low'
    benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load);
    
    toc_load = 'high'
    benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load);
    
    
    %% Middelburg SAR  -- only do best estimates here
    SAR='Middelburg';
    string_out = 'DOU_calc_Middelburg';
    %
    % toc_load = 'best'
    toc_load = 'low'
    res_Middelburg_best = benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load);
    toc_load = 'high'
    res_Middelburg_best = benthic_test.calc_a_from_Jorgensen_DOU(SAR, Db_Middelburg, string_out, 800, toc_load);
end