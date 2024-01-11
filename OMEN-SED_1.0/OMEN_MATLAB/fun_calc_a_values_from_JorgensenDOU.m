

% calls this one here: 
% [res] = calc_a_from_Jorgensen_DOU(SAR_Restreppo, Db_Middelburg, string_out, Zinf, toc_load)


Db_Middelburg = false;  % we use our new fit, see /Documents/Projects/18_RECCAP2/OMEN/Bioturbation/Solan_ea_2019/make_empirical_fct.m
Zinf = 800;     % enough fo r

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% orginal toc from MArkus without back calculation


%% Restreppo SAR
SAR_Restreppo = true;
string_out = 'DOU_calculation_Restrep';

% toc_load = 'best'
% res_Restreppo_best = benthic_test.calc_a_from_Jorgensen_DOU(SAR_Restreppo, Db_Middelburg, string_out, 800, toc_load);

 toc_load = 'low'
 res_Restreppo_low = benthic_test.calc_a_from_Jorgensen_DOU(SAR_Restreppo, Db_Middelburg, string_out, 800, toc_load);
 
 toc_load = 'high'
 res_Restreppo_high = benthic_test.calc_a_from_Jorgensen_DOU(SAR_Restreppo, Db_Middelburg, string_out, 800, toc_load);


%% Burwicz SAR
SAR_Restreppo=false;
string_out = 'DOU_calculation_Burwicz';

% toc_load = 'best'
% res_Burwicz_best = benthic_test.calc_a_from_Jorgensen_DOU(SAR_Restreppo, Db_Middelburg, string_out, 800, toc_load);

 toc_load = 'low'
 benthic_test.calc_a_from_Jorgensen_DOU(SAR_Restreppo, Db_Middelburg, string_out, 800, toc_load);
 
 toc_load = 'high'
 benthic_test.calc_a_from_Jorgensen_DOU(SAR_Restreppo, Db_Middelburg, string_out, 800, toc_load);


