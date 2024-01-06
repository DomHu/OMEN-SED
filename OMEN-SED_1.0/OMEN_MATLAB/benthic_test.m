%% OMEN-SED 1.1 BENTHIC-MODEL Stand-alone matlab code
% Hülse et al. (2017) GMD
% including the RCM approximation of Pika et al. (2020) GMD

% benthic_test.m
% functions to run OMEN-SED and plot the results

% Command to run the model:
% benthic_test.run_OMEN: For original 2G model
% benthic_test.run_OMEN_RCM: For RCM approximation
% 1) default sediment-water interface boundary conditions are set as prescribed in default_swi()
% 2) the subroutines for the different properties are called in test_benthic(1,swi):
% 3) results saved in res and are plotted with plot_column()

classdef benthic_test
    % test cases for benthic layer model
    
    properties
    end
    
    methods(Static)
        
        function swi = default_swi()
            % set default SWI conditions
            
            format longEng
            
            bsd = benthic_main();
            
            swi.Test_Dale = false;
            swi.Test_Dale_14G = false;               % use for 14G model as in Dale or nG as specified below
            swi.plot_fig = true;                                % plot the sediment profiles
            swi.write_output = true;
            
            swi.calc_P_DIC_ALK=false;       % also calculate P, DIC & ALK?
            
            swi.IntConst_GMD= false;         % true A1, A2 as in GMD; false: use Sandra's calculation
            
            swi.Nitrogen=true;                                  % calculate N (true/false)
            swi.Iron=false;                                      % calculate Fe (true/false)
            
            % for 2G-model
            swi.C01_nonbio= 0.7*1e-2/12*bsd.rho_sed;            % TOC concentration at SWI (wt%) -> (mol/cm^3 solid phase)
            swi.C02_nonbio= 0.3*1e-2/12*bsd.rho_sed;            % TOC concentration at SWI (wt%) -> (mol/cm^3 solid phase)
            swi.Fnonbio1 =  swi.C01_nonbio*(1-bsd.por)*bsd.w; %0.5*6.2e-3/100^2*365;    % swi.C01_nonbio*(1-bsd.por)*bsd.w;    % calculate flux [mol/(cm2 yr)] according non-bioturbated flux
            swi.Fnonbio2 =  swi.C02_nonbio*(1-bsd.por)*bsd.w; %0.5*6.2e-3/100^2*365;    % swi.C02_nonbio*(1-bsd.por)*bsd.w;    % calculate flux [mol/(cm2 yr)] according non-bioturbated flux
            swi.C01 = swi.C01_nonbio;                           % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
            swi.C02 = swi.C02_nonbio;                           % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
            
            % for nG-model
            swi.nG = 500;
            swi.p_a = 0.1;
            swi.p_nu = 0.151;
            swi.C0 = 1.0 * 1e-2/12*bsd.rho_sed;                 % TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
            swi.TOCwt_SWI = 4.5;
            swi.C0_nonbio = NaN; %swi.TOCwt_SWI * 1e-2/12*bsd.rho_sed;                 % TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)

            swi.FOM_total = NaN; %3.91140e-003;                                       % In case a total settling flux of OM is, specify this here and set swi.flux = true; [Note: not tested] 
            swi.flux = false;
            swi.k = 0.0;    % in case 1G as Thullner ea. '09
            
            %            swi.FeIII0=2.8E-005; %3.0E-006;                   	% FeIII concentration at SWI (mol/cm^3) --> TODO: needs to be a flux!
            Fe_influx_Dale = 0.1667*1110.0E-006;                     % FeIII influx from Dale but just 16.67% is available for DIR (See his Tab. 2, footnote d)
            swi.Flux_FeIII0 =  (Fe_influx_Dale)*365/100^2;           % Dale 1110 mol/(m^2 day)   -->  mol/(cm^2 yr):     *365/100^2
            %            swi.FeIII0= swi.Flux_FeIII0/((1-bsd.por)*bsd.w);     % calculate concentration [mol/cm^3] from flux [mol/(cm2 yr)] according non-bioturbated flux!!!
            
            %bottom water concentrations
            swi.T = 12.5;                                        % temperature (degree C)
            swi.O20=210.0E-009;                                 % O2  concentration at SWI (mol/cm^3)
            swi.NO30=9.6e-9;                                   % NO3 concentration at SWI (mol/cm^3)
            swi.NH40=0.4e-9;                                 	% NH4 concentration at SWI (mol/cm^3)
            swi.Fe20=0.0;                                       % Fe2 concentration at SWI (mol/cm^3)
            swi.SO40=2.8E-005;                                	% SO4 concentration at SWI (mol/cm^3)
            swi.H2S0=0.0;                                       % H2S concentration at SWI (mol/cm^3)
            swi.PO40=40.0e-9;                                   % PO4 concentration at SWI (mol/cm^3)
            swi.Mflux0=365*0.2e-10*1/(1-bsd.por)*1/bsd.w;       % actually CONCENTRATION of M at the sediment [mol/cm3] : from flux input  365*0.2e-10 (mol/(cm2*yr))
            swi.DIC0=2.4E-006;                                 	% DIC concentration at SWI (mol/cm^3)
            swi.ALK0=2.4E-006;                                 	% ALK concentration at SWI (mol/cm^3)
            swi.S0=35;                                         	% Salinity at SWI (not used at the moment)
            
            
            swi.plot_PO4_DIC_ALK=false;
        end
        
        
    	function [res] = calc_a_from_Jorgensen_DOU(SAR_Restreppo, Db_Middelburg, string_out)
            %% Calculate a values based on the DOU map of Jorgensen et al. (2022)
            %% Use all our updated boundary conditions, e.g., the Restreppo sedimentation rates, spatially explicit porosity values, zbio after Song et al., and new surface TOC wt%. 
            
            
            % SAR_Restreppo:    true : use Restreppo data 
            %                   false : use Burwicz parameterization
            % Db_Middelburg:    true : use Middelburg parametrization
            %                   false : use Dbio extrapolated from Solan data (not good) 
            
            zbio_fix = false;   % set this true to use zbio_global everywhere
            zbio_global=5.75;   % 5.75 cm (Teal et al., 2008) and 9.8 cm (Boudreau, 1998)
                        
            warning('query','all')
            warning('off','all')
            warning
            
            str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
            
            string_out = [string_out  'Best_estimate_'];
            SA_coefficient = 1.0;
            
            lat = (89.875:-0.25:-89.875);
            long = (-179.875:0.25:179.875);
 
            %__________________________________________________________________________
            %   load data
            %__________________________________________________________________________
            toc = load('./data/BC_calc_a_from_Jorgensen/surfOC_matrix_2023_new_v01.csv');
            porosity_matrix_new = load('./data/BC_calc_a_from_Jorgensen/porosity_matrix_2023_new_v01.csv');
            
            % load Restreppo SAR
            SAR_Restreppo_data = load('./data/BC_calc_a_from_Jorgensen/Restreppo_sar_matrix_2023-09-21_new.csv');            
            
            water_depth = load('./data/BC_calc_a_from_Jorgensen/bathymetry_matrix_new.csv');
            
            DOU_mask_coast_filled = struct2array(load('./data/BC_calc_a_from_Jorgensen/DOU_mask_coast_filled.mat'));
            
            % zbio interpolated from Song et al. (2022)
            zbio_interpolated = struct2array(load('./data/BC_calc_a_from_Jorgensen/Bioturbation_depth_remapped_natural.mat'));
            zbio_nan=0;
            zbio_good=0;
            
            % load a-values from Pika et al. (2023) to start with
         	a_values_start = struct2array(load('./data/BC_calc_a_from_Jorgensen/a_value_local_Start.mat'));

            %% TODO FROM HERE - load the updated BCs from ./data/BC_calc_a_from_Jorgensen  -- check with run_OMEN_RCM_global how things are called!
            BC.Tmp = struct2array(load('./data/BC_calc_a_from_Jorgensen/BC_Tmp_Margin_updated.mat'));                             % BW temperature [°C]
            BC.O2 = struct2array(load('./data/BC_calc_a_from_Jorgensen/BC_O2_Margin_updated.mat'));                   	% BW  O2 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
            BC.NO3 = struct2array(load('./data/BC_calc_a_from_Jorgensen/BC_NO3_Margin_updated.mat'));                             % BW  NO3 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
            

            % error counters
            wrong_iter = 0;
          	wrong_iter_age = 0;
            zbioMatching_xy_GMD(1,:)=[0 0];
            zbioMatching_xy_final(1,:)=[0 0];
         	NEGATIVE_OXID_xy(1,:)=[0 0];

            debug.iter_zbioMatching_GMD = 0;        % counter for calculating the correct TOC profile using the GMD IntConst
           	debug.iter_zbioMatching_final = 0;      % counter for calculating the correct TOC profile after trying different methods
           	debug.iter_NEGATIVE_OXID = 0;           % counter for negative total OM oxidation


            [m,n]=size(toc);

            ncl = 1;
            res.bsd = benthic_main(ncl);            
         	res.bsd.zinf = 800;	

            
            res.bsd.usescalarcode = ncl==1;
            
            swi = benthic_test.default_swi();
            swi.TwoG_OM_model = false;
         	swi.flux = false;

            
            %__________________________________________________________________________
            %   set specific parameters
            %__________________________________________________________________________
            
            res.swi = swi;
            
            xstart=1; %55;
            xstop=m;  %55;        %lat
            ystart=1;  % 6;
            ystop=n;  % 6;        %long (or the other way around- who knows!)

            Restreppo_NaN = 0;
            for x = 1:m
%         	for x = xstart:xstop
                x
                % calculate volume
                % convert deg to cm CODAS package (by E.Firing,et al.)
                rlat = lat(x) * pi/180;
                m = 111132.09  - 566.05 * cos(2 * rlat)+ 1.2 * cos(4 * rlat);
                dy = 0.25*m*100.0; %cm
                p = 111415.13 * cos(rlat) - 94.55 * cos(3 * rlat);
                dx = 0.25*p*100.0; %cm
                
                for y = 1:n
%                for y = ystart:ystop
                    if ((isnan(toc(x,y))))    % check for toc = NaN

                      	dxdy(x,y)   = NaN;                                               
                        a_value_best(x,y) = NaN;

                    else
                        
                        % Here normal call: RUN OMEN-SED 
                        % res=benthic_test.test_benthic(1,swi);
                        % TODO: Do a new function (similar to benthic_test.test_benthic) 
                        % for fitting DOU values with all the below in it
                        % First give it all swi boundary consitions, see run_OMEN_RCM_global(exp_name, Zinf)
                        
                        % in RECCAP2:                        
                        res.swi.p_a=a_values_start(x,y);    % set first guess for a-value
                        res.bsd.wdepth = -water_depth(x,y);
                        if(SAR_Restreppo)
                            res.bsd.w=SAR_Restreppo_data(x,y)*SA_coefficient;
                            if(isnan(res.bsd.w))
                                res.bsd.w=benthic_main.sedrate(res.bsd.wdepth)*SA_coefficient;
                                Restreppo_NaN = Restreppo_NaN+1;
                            end
                        else
                            res.bsd.w=benthic_main.sedrate(res.bsd.wdepth)*SA_coefficient;
                        end                        

                        if(Db_Middelburg)
                            res.bsd.Dbio=benthic_main.biorate(res.bsd.wdepth)*SA_coefficient;                       
                        else
                            res.bsd.Dbio=benthic_main.biorate_Solan_MARCATS(res.bsd.wdepth)*SA_coefficient;    
                        end
                        res.bsd.zbio=zbio_interpolated(x,y)*SA_coefficient;
                        
                        if(zbio_fix)
                            res.bsd.zbio=zbio_global;   % 5.75 cm (Teal et al., 2008) and 9.8 cm (Boudreau, 1998)
                        end
                        if(isnan( res.bsd.zbio))
                             res.bsd.zbio=zbio_global*SA_coefficient;  % set to zbio_global if no value
                             zbio_nan=zbio_nan+1;
                        else
                            zbio_good=zbio_good+1;
                        end
                        
                        res.bsd.por=porosity_matrix_new(x,y)/100;  %porosity at SWI
                        
                        % initialize parameters of analytical solution
%                        rho_sed_loc = 2.5;
                        conv=rho_sed/100*(1-res.bsd.por);   % convert wt% -> g/cm3 (total sediment)   -- this is what Sandra did for James
                        conv=rho_sed/100;   % convert wt% -> g/cm3 (total sediment)     -- I think in my model I need it as bulk sediment because I apply *(1-por) when I calculate the flux in calcCflx() !!??
                        conv=res.bsd.rho_sed/(100*12);   % convert wt% -> mol/cm3 (total sediment)     -- I think in my model I need it as bulk sediment because I apply *(1-por) when I calculate the flux in calcCflx() !!??
                        res.swi.C0=toc(x,y)*conv;          % POC at SWI (wt% -> mol/cm3 (total sediment))
                        
                        
                         % initialize & calculate
                        res.zTOC_RCM = benthic_zTOC_RCM(res.bsd);
                        % Adding into on RCM for MultiG approach
                        [res.zTOC_RCM.k, res.swi.C0i, res.swi.Fnonbioi] = benthic_test.RCM(res.bsd, res.swi);
                        res = res.zTOC_RCM.calc(res.bsd,res.swi, res);
                        
                        dxdy(x,y)   = dx*dy;    % cm^2 per grid-cell

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% calculate TOC fluxes 
                     	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %%%%%%%%%%%%%%%%%
                        % flux through SWI (just for BE calculation - total value is way too high -- maybe a few very wrong grid-cells!?)
                        [F_TOC_swi, F_TOC1_swi] = res.zTOC_RCM.calcCflx(0, res.bsd, res.swi, res);
                        Flux_TOC_swi(x,y) = nansum(F_TOC1_swi)*dx*dy;  % 
                        Flux_TOC_swi_perArea(x,y) = nansum(F_TOC1_swi);  % 

                        
                        %%%%%%%%%%%%%%%%%
                        % through different depth levels - calculate maps of BE and Fluxes                         
                        for i=1:length(depth_levels)

                            [F_TOC_inf, F_TOC1_inf] = res.zTOC_RCM.calcCflx(depth_levels(i), res.bsd, res.swi, res);
                            if(nansum(F_TOC1_inf)<nansum(F_TOC1_swi) | nansum(F_TOC1_inf) >0 )   % check if flu in sedimentis larger at SWI
                              %  nansum(F_TOC1_inf);
                                wrong_iter=wrong_iter+1;
                                F_TOC1_inf=0;
                            end
                            Flux_depth{i}(x,y) = nansum(F_TOC1_inf)*dx*dy;  %   in g/yr
                         	Flux_depth_inTg{i}(x,y) = Flux_depth{i}(x,y)*1e-12;   %   in Tg/yr
                            Flux_depth_perArea{i}(x,y) = nansum(F_TOC1_inf);  %   in g /(cm^2 yr)
                            BE_depth{i}(x,y) = nansum(F_TOC1_inf)/nansum(F_TOC1_swi)*100;                            
                        end

                        % through different age-levels - calculate maps of BE and Fluxes             
                        for j=1:length(age_levels)
                            %  age_levels = [150 1000 5000];  % in yr    
                            % NOTE: add zbio as the age related depth-layers start with 0 = zbio
                            [F_TOC_inf, F_TOC1_inf] = res.zTOC_RCM.calcCflx(Age_layers{j}(x,y) +  res.bsd.zbio, res.bsd, res.swi, res);
                            if(nansum(F_TOC1_inf)<nansum(F_TOC1_swi) | nansum(F_TOC1_inf) > 0 )   % check if flu in sedimentis larger at SWI
                              %  nansum(F_TOC1_inf);
                                wrong_iter_age=wrong_iter_age+1;
                                F_TOC1_inf=0;
                            end
                            Flux_age{j}(x,y) = nansum(F_TOC1_inf)*dx*dy;  %   in g/yr
                         	Flux_age_inTg{j}(x,y) = Flux_age{j}(x,y)*1e-12;   %   in Tg/yr
                            Flux_age_perArea{j}(x,y) = nansum(F_TOC1_inf);  %   in g /(cm^2 yr)
                            BE_age{j}(x,y) = nansum(F_TOC1_inf)/nansum(F_TOC1_swi)*100;                            
                        end                        

                    end
                    
                    
                end
            end
            wrong_iter
            wrong_iter_age
           	Restreppo_NaN

            % calculate global total burial per depth level  in Tg/yr
            for i=1:length(depth_levels)
                Total_burial_depth(i) = nansum(nansum(Flux_depth{i}))/10^(12);
            end
            % calculate global total burial per age level  in g/yr
            for i=1:length(age_levels)
                Total_burial_age(i) = nansum(nansum(Flux_age{i}))/10^(12);
            end                    
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % save results
            save(['output/Flux_TOC_swi_', string_out, str_date '.mat'] , 'Flux_TOC_swi');            

%            save(['output/BE_depth_', string_out, str_date '.mat'] , 'BE_depth');
            save(['output/' str_date '_Flux_depth_', string_out '.mat'] , 'Flux_depth');                   %   in g/yr
            save(['output/' str_date '_Flux_depth_inTg_', string_out '.mat'] , 'Flux_depth_inTg');         %   in Tg/yr
            save(['output/' str_date '_Flux_depth_perArea', string_out '.mat'] , 'Flux_depth_perArea');    %   in g /(cm^2 yr)                      
         	save(['output/' str_date '_Total_burial_depth_', string_out '.mat'] , 'Total_burial_depth');
            
%         	save(['output/BE_age_', string_out, str_date '.mat'] , 'BE_age');
            save(['output/' str_date '_Flux_age_', string_out '.mat'] , 'Flux_age');                   %   in g/yr
            save(['output/' str_date '_Flux_age_inTg_', string_out '.mat'] , 'Flux_age_inTg');         %   in Tg/yr
            save(['output/' str_date '_Flux_age_perArea', string_out '.mat'] , 'Flux_age_perArea');    %   in g /(cm^2 yr)                      
         	save(['output/' str_date '_Total_burial_age_', string_out '.mat'] , 'Total_burial_age');
            
            %%  Calculate the TOC burial fluxes for different areas  (i.e., MARCAT, EEZ, Longhurst)            
            
            % 45 MARCAT regions (only ~39 with TOC values)
            load('MARCAT_Regions_UPDATE/TOC_polygon_MARCAT.mat');
          	S = shaperead('./MARCAT_Regions_UPDATE/Regions_V2.shp');
            for i=1:length(S)
                Mean_lat(i)=nanmean(S(i).Y);
            end

            if true % was  false
            [TOC_burial_per_OceanRegion, Area_per_OceanRegion] = benthic_test.calc_area_burial(Flux_depth_inTg, dxdy, TOC_polygon_MARCAT);
            Tr_TOC_burial_per_OceanRegion= TOC_burial_per_OceanRegion';
            Tr_Area_per_OceanRegion=Area_per_OceanRegion';
            Table_Ocean_region = array2table([vertcat(S.OBJECTID) Tr_TOC_burial_per_OceanRegion Tr_Area_per_OceanRegion, Tr_TOC_burial_per_OceanRegion(:,5)./Tr_Area_per_OceanRegion(:).*10^9, Mean_lat', vertcat(S.MARCAT)]);
            Table_Ocean_region.Properties.VariableNames(1:12) = {'OBJECTID','1cm','10cm','20cm','50cm', '100cm', '200cm','500cm', 'Size (km2)', 'Area flx 1m','Mean lat','MARCAT'};
            save(['output/' str_date '_Table_OceanRegion_MARCAT_45regions', string_out '.mat'] , 'Table_Ocean_region');
            end

            [TOC_burial_per_OceanRegion_depth, Area_per_OceanRegion_depth] = benthic_test.calc_area_burial(Flux_age_inTg, dxdy, TOC_polygon_MARCAT);
            Tr_TOC_burial_per_OceanRegion_depth= TOC_burial_per_OceanRegion_depth';
            Tr_Area_per_OceanRegion_depth=Area_per_OceanRegion_depth';
            Table_Ocean_region_age = array2table([vertcat(S.OBJECTID) Tr_TOC_burial_per_OceanRegion_depth Tr_Area_per_OceanRegion_depth, Mean_lat', vertcat(S.MARCAT)]);
            Table_Ocean_region_age.Properties.VariableNames(1:7) = {'OBJECTID','100years','500years','1000years', 'Size (km2)','Mean lat','MARCAT'};
            save(['output/' str_date '_Table_OceanRegion_MARCAT_45regions_age', string_out '.mat'] , 'Table_Ocean_region_age');
            
            if true % was false
            % MARCAT areas:
            % 1: EBC, 2: WBC, 3: Indian, 4: Sub-Polar, 5: Polar, 6: Marginal, 7: Tropical                        
            load('MARCAT_Regions_UPDATE/TOC_MARCAT.mat', 'TOC_MARCAT');
            
            if  true % was false
            [TOC_burial_per_OceanRegion, Area_per_OceanRegion] = benthic_test.calc_area_burial(Flux_depth_inTg, dxdy, TOC_MARCAT);
            Tr_TOC_burial_per_OceanRegion= TOC_burial_per_OceanRegion';
            Tr_Area_per_OceanRegion=Area_per_OceanRegion';
            Table_Ocean_region = array2table([(1:7)' ,Tr_TOC_burial_per_OceanRegion Tr_Area_per_OceanRegion, Tr_TOC_burial_per_OceanRegion(:,5)./Tr_Area_per_OceanRegion(:).*10^9]);
            Table_Ocean_region.Properties.VariableNames(1:10) = {'MARCAT','1cm','10cm','20cm','50cm', '100cm', '200cm','500cm', 'Size (km2)', 'Area flx 1m'};
            save(['output/' str_date '_Table_OceanRegion_MARCAT', string_out '.mat'] , 'Table_Ocean_region');
            end
            [TOC_burial_per_OceanRegion_age, Area_per_OceanRegion_age] = benthic_test.calc_area_burial(Flux_age_inTg, dxdy, TOC_MARCAT);
            Tr_TOC_burial_per_OceanRegion_age= TOC_burial_per_OceanRegion_age';
            Tr_Area_per_OceanRegion_age=Area_per_OceanRegion_age';
            Table_Ocean_region_age = array2table([(1:7)' ,Tr_TOC_burial_per_OceanRegion_age Tr_Area_per_OceanRegion_age]);
            Table_Ocean_region_age.Properties.VariableNames(1:5) = {'MARCAT','100years','500years','1000years', 'Size (km2)'};
            save(['output/' str_date '_Table_OceanRegion_MARCAT_age', string_out '.mat'] , 'Table_Ocean_region_age');
            
            
            % Longhurst areas:
          	load('longhurst_v4_2010/TOC_polygon_Longhurst.mat')
            S = shaperead('./longhurst_v4_2010/Longhurst_world_v4_2010.shp');
 
            [TOC_burial_per_OceanRegion, Area_per_OceanRegion] = benthic_test.calc_area_burial(Flux_depth_inTg, dxdy, TOC_polygon_Longhurst);            
            Tr_TOC_burial_per_OceanRegion= TOC_burial_per_OceanRegion';
            Tr_Area_per_OceanRegion=Area_per_OceanRegion';
           	Table_Ocean_region = array2table([(1:length(S))' Tr_TOC_burial_per_OceanRegion Tr_Area_per_OceanRegion, Tr_TOC_burial_per_OceanRegion(:,5)./Tr_Area_per_OceanRegion(:).*10^9]);
            Table_Ocean_region.Properties.VariableNames(1:10) = {'OBJECTID','1cm','10cm','20cm','50cm', '100cm', '200cm','500cm', 'Size (km2)', 'Area flx 1m'};
            save(['output/' str_date '_Table_OceanRegion_Longhurst', string_out '.mat'] , 'Table_Ocean_region');

            [TOC_burial_per_OceanRegion_age, Area_per_OceanRegion_age] = benthic_test.calc_area_burial(Flux_age_inTg, dxdy, TOC_polygon_Longhurst);            
            Tr_TOC_burial_per_OceanRegion_age= TOC_burial_per_OceanRegion_age';
            Tr_Area_per_OceanRegion_age=Area_per_OceanRegion_age';
           	Table_Ocean_region_age = array2table([(1:length(S))' Tr_TOC_burial_per_OceanRegion_age Tr_Area_per_OceanRegion_age]);
            Table_Ocean_region_age.Properties.VariableNames(1:5) = {'OBJECTID','100years','500years','1000years', 'Size (km2)'};
            save(['output/' str_date '_Table_OceanRegion_Longhurst_age', string_out '.mat'] , 'Table_Ocean_region_age');
            

            % EEZ
            S_EEZ = shaperead('./EEZ/eez_v11.shp');
            
            % EEZ countries code 1:157
          	load('EEZ/TOC_EEZ_Country.mat')
            
            [TOC_burial_per_OceanRegion, Area_per_OceanRegion] = benthic_test.calc_area_burial(Flux_depth_inTg, dxdy, TOC_EEZ_Country);
            Tr_TOC_burial_per_OceanRegion= TOC_burial_per_OceanRegion';
            Tr_Area_per_OceanRegion=Area_per_OceanRegion';            
           	Table_Ocean_region = array2table([(1:157)' Tr_TOC_burial_per_OceanRegion Tr_Area_per_OceanRegion, Tr_TOC_burial_per_OceanRegion(:,5)./Tr_Area_per_OceanRegion(:).*10^9]);
            Table_Ocean_region.Properties.VariableNames(1:10) = {'Country_ID','1cm','10cm','20cm','50cm', '100cm', '200cm','500cm', 'Size (km2)', 'Area flx 1m'};
            save(['output/' str_date '_Table_OceanRegion_EEZ', string_out '.mat'] , 'Table_Ocean_region');
  
            [TOC_burial_per_OceanRegion_age, Area_per_OceanRegion_age] = benthic_test.calc_area_burial(Flux_age_inTg, dxdy, TOC_EEZ_Country);
            Tr_TOC_burial_per_OceanRegion_age= TOC_burial_per_OceanRegion_age';
            Tr_Area_per_OceanRegion_age=Area_per_OceanRegion_age';            
           	Table_Ocean_region_age = array2table([(1:157)' Tr_TOC_burial_per_OceanRegion_age Tr_Area_per_OceanRegion_age]);
            Table_Ocean_region_age.Properties.VariableNames(1:5) = {'Country_ID','100years','500years','1000years', 'Size (km2)'};
            save(['output/' str_date '_Table_OceanRegion_EEZ_age', string_out '.mat'] , 'Table_Ocean_region_age');
  
            end
            % This is all tiny EEZ territories (above they are combined per country)
            if false    
                % EEZ territories:
                load('EEZ/TOC_polygon_EEZ.mat')

                [TOC_burial_per_OceanRegion, Area_per_OceanRegion] = benthic_test.calc_area_burial(Flux_depth_inTg, dxdy, TOC_polygon_EEZ);
                Tr_TOC_burial_per_OceanRegion= TOC_burial_per_OceanRegion';
                Tr_Area_per_OceanRegion=Area_per_OceanRegion';
                Table_Ocean_region = array2table([(1:length(S_EEZ))' Tr_TOC_burial_per_OceanRegion Tr_Area_per_OceanRegion, Tr_TOC_burial_per_OceanRegion(:,5)./Tr_Area_per_OceanRegion(:).*10^9]);
                Table_Ocean_region.Properties.VariableNames(1:10) = {'TERRITORY_ID','1cm','10cm','20cm','50cm', '100cm', '200cm','500cm', 'Size (km2)', 'Area flx 1m'};
                save(['output/' str_date '_Table_OceanRegion_EEZ_Territory', string_out '.mat'] , 'Table_Ocean_region');
            end
            
            zbio_nan
            zbio_good
        end

        
        
        function run_OMEN()
            % run OMEN-SED with default SWI conditions as in default_swi()
            
            tic;
            swi=benthic_test.default_swi();
            swi.TwoG_OM_model = true;
                        swi.BC_wdepth_flag = false;
            swi.BC_sed_rate_flag = false;
            swi.flux = false;
            swi.plot_fig = true;                                % plot the sediment profiles
%            % set date-time
            %            str_date = datestr(now,'ddmmyy_HH_MM_SS');
            res=benthic_test.test_benthic(1,swi);
            toc;
            benthic_test.plot_column(res, false, res.swi, '0706_2G_OMEN')
            
            % calculate depth integrated OM degradation rates (wrong use
            % difference of fluxes (see test_benthic() below)
            Cox_rate.Cox_total = res.zTOC.calcReac(0.0, res.bsd.zinf, 1*(1-res.bsd.por), 1*(1-res.bsd.por), res.bsd, res.swi, res);
            Cox_rate.Cox_Aerobic = res.zTOC.calcReac(0.0, res.zox, 1*(1-res.bsd.por), 1*(1-res.bsd.por), res.bsd, res.swi, res);
            if(res.swi.Nitrogen)
                Cox_rate.Cox_Denitr = res.zTOC.calcReac(res.zox, res.zno3, 1*(1-res.bsd.por), 1*(1-res.bsd.por), res.bsd, res.swi, res);
            end
            if(res.swi.Iron)
                Cox_rate.Cox_IronIII = res.zTOC.calcReac(res.zno3, res.zfeIII, 1*(1-res.bsd.por), 1*(1-res.bsd.por), res.bsd, res.swi, res);
            end
            Cox_rate.Cox_sulfate = res.zTOC.calcReac(res.zno3, res.bsd.zinf, 1*(1-res.bsd.por), 1*(1-res.bsd.por), res.bsd, res.swi, res)
            
            % Calculate out put to compare with Dale:
            %            Cox_total_Dale_units = Cox_rate.Cox_total*1000 *100^2/365      % in units of mmol m-2 d-1 as in Dale ea. 2015, Fig. 2a
            
            % calculate mean OM concentration in upper x cm
            [C_10, C1_10, C2_10] = res.zTOC.calcC( 10, res.bsd, res.swi, res);
            OM_10=C_10* 100*12/res.bsd.rho_sed;
            x = 10;
            Mean_OM = 1/x * 100*12/res.bsd.rho_sed*res.zTOC.calcOM(0.0, x, 1, 1, res.bsd, res.swi, res);
        end
        
        function [Cox_rate, res, Output]  = run_OMEN_RCM(exp_name)
            % run OMEN-SED with default SWI conditions as in default_swi()
            
            swi = benthic_test.default_swi();
            swi.TwoG_OM_model = false;
            swi.BC_wdepth_flag = false;
            swi.BC_sed_rate_flag = false;

            swi.plot_fig = true;                                % plot the sediment profiles
            
            if(swi.Test_Dale)
                if(swi.Test_Dale_14G)
                    swi.nG = 14;
                end
                swi.p_a = 0.1;
                swi.p_nu = 0.125;
                swi.POC_flux = [0.5 1 2 4 6 8 10 12 14 16];
                swi.POCi = 1;
                O20_v = [1e-9 2e-9 5e-9 10e-9 15e-9 25e-9 50e-9 100e-9 200e-9];
                swi.O20=O20_v(3);
            end
            
            res=benthic_test.test_benthic(1,swi);
            % set date-time or string going into plot function
            str_date = [exp_name, num2str(res.swi.nG) 'G_a=' num2str(res.swi.p_a) '_nu=' num2str(res.swi.p_nu) '_O20_' num2str(res.swi.O20)];
            if(swi.plot_fig)
                benthic_test.plot_column(res, false, swi, str_date)
                benthic_test.plot_TOC(res, false, swi, str_date)
            end
            % calculate depth integrated OM degradation rates (this is
            % wrong - use the fluxes of TOC at the boundaries - see test_benthic)
            swi.C0i = res.swi.C0i;
            Cox_rate(1,1) = res.zTOC_RCM.calcReac(0.0, res.bsd.zinf, 1*(1-res.bsd.por), res.bsd, swi, res);
            Cox_rate(1,2) = res.zTOC_RCM.calcReac(0.0, res.zox, 1*(1-res.bsd.por), res.bsd, swi, res)/Cox_rate(1,1)*100;
            if(swi.Nitrogen)
                Cox_rate(1,3) = res.zTOC_RCM.calcReac(res.zox, res.zno3, 1*(1-res.bsd.por), res.bsd, swi, res)/Cox_rate(1,1)*100;
            end
            if(res.swi.Iron)
                Cox_rate(1,4) = res.zTOC_RCM.calcReac(res.zno3, res.zfeIII,1*(1-res.bsd.por), res.bsd, swi, res)/Cox_rate(1,1)*100;
            end
            Cox_rate(1,5) = res.zTOC_RCM.calcReac(res.zfeIII,res.zso4, 1*(1-res.bsd.por), res.bsd, swi, res)/Cox_rate(1,1)*100;
            
            if(swi.Test_Dale)
                % compare results with Dale
                [F_TOC_swi, F_TOC1_swi] = res.zTOC_RCM.calcCflx(0, res.bsd, res.swi, res);
                [F_TOC_inf, F_TOC1_inf] = res.zTOC_RCM.calcCflx(res.bsd.zinf, res.bsd, res.swi, res);
                
                Cox_total = (F_TOC_swi-F_TOC_inf)*1000 *100^2/365;      % in units of mmol m-2 d-1 as in Dale ea. 2015, Fig. 2a
                Cox_rate.Flux_Fe2_Dale_units = res.flxswiFe2 *10^6 *100^2 /365;                   % in units of umol m-2 d-1 as in Dale ea. 2015, Fig. 3
                fprintf('\n');
                fprintf('Total Cox as flux difference (mmol m-2 d-1) %e \n',  Cox_total);
                fprintf('SWI-flux FeII (umol m-2 d-1) %e \n',  Flux_Fe2_Dale_units);
                
            end
            
            
            % calculate mean OM concentration in upper x cm
            [C_10, C1_11] = res.zTOC_RCM.calcC( 10, res.bsd, res.swi, res);
            OM_10=C_10* 100*12/res.bsd.rho_sed;
            x = 10;
            Mean_OM = 1/x * 100*12/res.bsd.rho_sed*res.zTOC_RCM.calcOM(0.0, x, 1, res.bsd, swi, res);
            
            % SWI-fluxes in mol cm-2 yr-1
            Output(1,1) = res.flxswiH2S;
            Output(1,2) = res.flxswiSO4;
            if(swi.Iron)
                Output(1,3) = res.flxswiFe2;
                Output(1,4)  = res.flxswiFe2 *10^6 *100^2 /365;                   % in units of umol m-2 d-1 as in Dale ea. 2015, Fig. 3
                Output(1,5) = res.flxswiFeIII;
            end
            Output(1,6) = res.flxswiNH4;
            Output(1,7) = res.flxswiNO3;
            Output(1,8) = res.flxswiO2;
            Output(1,9) = res.zox;
            Output(1,10) = res.zno3;
            Output(1,11) = res.zfeIII;
            Output(1,12) = res.zso4;
        end
        
        function [TOC_SWI_BC, TOC_MEAN_5cm] = inversely_calc_TOC_BC_global(exp_name, Zinf)
            %% inversely calculate the SWI BC for TOC so that the mean wt% in upper 5cm is the same as Lee et al. (2019)
            
            calc_res = 2;           % 1: 1/4°;  2: 1°;  3: 2°
            
            swi = benthic_test.default_swi();
            swi.flux = false;
            
            swi.sed_column_depth = Zinf;  % use Zinf as bsd.zinf
            
            swi.TwoG_OM_model = false;
            OneGmodel = false;
            swi.Test_Dale = false;  % he used the same a-value evrywhere - I don't want that -- plus I use 100G
            
            switch calc_res
                case 1
                    
                    % to calculate volume in grid-cell
                    dlat = 0.25;    % latitude difference in degrees for high-res
                    dlon = 0.25;    % longitude difference in degrees for high-res
                    
                case 2
                    % load boundary conditions - all in 1 degree resolution
                    load('./data/01_plot_BC_April_21/BC_1degree/Lee_toc_lr_weighted.mat');       % TOC at SWI [wt%]     -- need to translate to mol/cm^3 solid phase, i.e. *1e-2/12*bsd.rho_sed
                    toc_Lee = Lee_toc_lr_weighted;
                    
                    %load('./data/01_plot_BC_April_21/BC_1degree/a_lr_aligned_1to100.mat');                                   %  Parameter a [yrs]  (the apparent initial age of the initial organic matter mixture ) as fct of sedimentation rate dependent formulation for a (Arndt et al., 2013)
                    %                load('./data/01_plot_BC_April_21/BC_1degree/a_lr_aligned_1to50.mat');                                   %  Parameter a [yrs]  (the apparent initial age of the initial organic matter mixture ) as fct of sedimentation rate dependent formulation for a (Arndt et al., 2013)
                    % a_hr_updated = a_lr_aligned;
                    load('./data/01_plot_BC_April_21/BC_1degree/a_lr_aligned_PK_filled_in_a1to10.mat');                                   %  Parameter a [yrs]  (the apparent initial age of the initial organic matter mixture ) as fct of sedimentation rate dependent formulation for a (Arndt et al., 2013)
                    a_hr_updated = Parameter_a_total_PKsedrate1to10;
                    load('./data/01_plot_BC_April_21/BC_1degree/sed_lr_weighted_aligned.mat');                           % Sedimentation rate [cm/yr] -- as in Bradley ea. 2020, after Burwicz ea. 2011
                    sed_holo_hr_updated = sed_lr_weighted_aligned;
                    load('./data/01_plot_BC_April_21/BC_1degree/water_depth_aligned.mat');                             % Seafloor depth [m]
                    water_depth_updated = -water_depth_aligned;
                    load('./data/01_plot_BC_April_21/BC_1degree/lat_lr.mat');
                    lat = lat_lr;
                    
                    % to calculate volume in grid-cell
                    dlat = 1.0;     % latitude difference in degrees for 1 degree resolution
                    dlon = 1.0;     % longitude difference in degrees for 1 degree resolution
                    
                case 3
                    % load boundary conditions - all in 2 degree resolution
                    % to calculate volume in grid-cell
                    dlat = 2.0;     % latitude difference in degrees for 2 degree resolution
                    dlon = 2.0;     % longitude difference in degrees for 2 degree resolution
                    
            end
            
            [m,n]=size(toc_Lee);
            tic
            %% loop through lat - lon
            xstart = 1;
            ystart = 1;
            for x=xstart:m
                fprintf('-------------- New latitude ----------- \n');
                fprintf('x, %i \n',  x);
                
                for y=ystart:n
                    %                 fprintf('\n');
                    fprintf('x, y, %i %i \n',  x, y);
                    % now check for sed_holo or toc = NaN -- if yes set output to NaN
                    if ((isnan(sed_holo_hr_updated(x,y))) | (isnan(toc_Lee(x,y))) | (isnan(a_hr_updated(x,y)))| (isnan(water_depth_updated(x,y))))
                        
                        TOC_SWI_BC(x,y) = NaN;
                        TOC_MEAN_5cm(x,y) = NaN;
                    else        % valid BC: run the model
                        
                        
                        rho_sed_loc = 2.5;
                        swi.C0_nonbio = toc_Lee(x,y)*1e-2/12*rho_sed_loc;
                        swi.p_a = a_hr_updated(x,y);
                        
                        swi.BC_sed_rate_flag = true;    % flag for existing BC for sed-rate
                        swi.BC_sed_rate = sed_holo_hr_updated(x,y);     % holocene sed-rate
                        
                        swi.BC_wdepth_flag = true;
                        swi.BC_wdepth = water_depth_updated(x,y);
                        if(swi.BC_wdepth>7900)
                            swi.BC_wdepth %=7900;
                        end
                        
                        res=benthic_test.calculate_meanTOC(1,swi);
                        % check that calculated mean TOC in upper 5cm is
                        % within 1% of value reported in Lee et al. (saved in toc_Lee)
                        toc_threshold = 0.01;
                        toc_SWI_BC_help = toc_Lee(x,y);
                        result_accurate = res.Mean_OM/toc_Lee(x,y);
                        
                        while( abs(result_accurate - 1.0)> toc_threshold )
                            % increase or decrease SWI TOC et%
                            toc_SWI_BC_help = toc_SWI_BC_help * (2.0 - res.Mean_OM/toc_Lee(x,y));   % cases below not needed, for both cases the same calculation is needed
                            swi.C0_nonbio = toc_SWI_BC_help*1e-2/12*rho_sed_loc;
                            res=benthic_test.calculate_meanTOC(1,swi);
                            result_accurate = res.Mean_OM/toc_Lee(x,y);
                            
                            %                            if(res.Mean_OM/toc_Lee(x,y) < 1.0)    % mean toc too small --> increase SWI-BC
                            %                                toc_SWI_BC_help = toc_SWI_BC_help * (2.0 - res.Mean_OM/toc_Lee(x,y));
                            %                            else        % mean toc too large --> decrease SWI-BC
                            %                                toc_SWI_BC_help = toc_SWI_BC_help * (2.0 - res.Mean_OM/toc_Lee(x,y));
                            %                            end
                        end
                        
                        TOC_SWI_BC(x,y) = toc_SWI_BC_help;
                        TOC_MEAN_5cm(x,y) = res.Mean_OM;
                        
                        
                    end
                end
            end
            a_min_new = nanmin(nanmin(a_hr_updated));
            a_max_new = nanmax(nanmax(a_hr_updated));
%             save(['./data/01_plot_BC_April_21/BC_1degree/Lee_toc_lr_weighted_SWI_a', num2str(a_min_new) , 'to', num2str(a_max_new) , exp_name ,'.mat'] , 'TOC_SWI_BC')
%             save(['./data/01_plot_BC_April_21/BC_1degree/Lee_toc_lr_weighted_mean5cm_a', num2str(a_min_new) , 'to', num2str(a_max_new) , exp_name ,'.mat'] , 'TOC_MEAN_5cm')
            save(['./data/01_plot_BC_April_21/BC_1degree/Lee_toc_lr_weighted_SWI_PK_filled_a1to10' , exp_name ,'.mat'] , 'TOC_SWI_BC')
            save(['./data/01_plot_BC_April_21/BC_1degree/Lee_toc_lr_weighted_mean5cm_PK_filled_a1to10' , exp_name ,'.mat'] , 'TOC_MEAN_5cm')
            
            
        end
        
        function [Cox_rate_out, Penetration_out, Flux_Fe2_Dale_units, debug, dxdy, SWI_fluxes, TOC_burial_flux] = run_OMEN_RCM_global(exp_name, Zinf)
            %% run global OMEN-SED with boundary conditions as in data
            % exp_name: string for experiment name
            % default_01: 1-degree resolution, 1m sediment-column, a-value = Arndt-parametr.
            
            
            warning('query','all')
            warning('off','all')
            warning
            
            depth_dep_a = false;     % use depth-dependent a-values?
            calc_res = 2;           % 1: 1/4°;  2: 1°;  3: 2°
            
            Fe3_HR_frac = 0.1667;           % 16.67% of total FeOOH-influx is HR and thus available for DIR (See his Tab. 2, footnote d)
            hypoxic_th = 63.0e-9;           % threshold for hypoxia after Dale et al. '15
            Fe3_in_SMI_frac_oxic = 0.6;     % Fe3_HR used for Sulphide-Mediated Iron reduction under oxic BW (i.e. not available for DIR)
            Fe3_in_SMI_frac_hypoxic = 0.3;  % Fe3_HR used for Sulphide-Mediated Iron reduction under hypoxic BW (i.e. not available for DIR)
            Cox_rate_out.grid_hypoxic = 0;
            
            Fe_fail_xy(1,:)=[0 0];
            zbioMatching_xy_GMD(1,:)=[0 0];
            zbioMatching_xy_final(1,:)=[0 0];
         	NEGATIVE_OXID_xy(1,:)=[0 0];

            
            swi = benthic_test.default_swi();
            swi.flux = false;
            
            swi.sed_column_depth = Zinf;  % use Zinf as bsd.zinf
            
            swi.TwoG_OM_model = false;
            OneGmodel = false;
            swi.Test_Dale = false;  % he used the same a-value evrywhere - I don't want that -- plus I use 100G
            %             swi.Test_Dale_14G = true;
            %            	if(swi.Test_Dale)
            %                 if(swi.Test_Dale_14G)
            %                     swi.nG = 14;
            %                 end
            %             end
            switch calc_res
                case 1
                    %            if(calc_high_res)
                    % load boundary conditions - all in 1/4 degree resolution
                    load './data/toc_in_NaN.mat'                                       % TOC at SWI [wt%]     -- need to translate to mol/cm^3 solid phase, i.e. *1e-2/12*bsd.rho_sed
                    load './data/a_hr_updated.mat'                                  %  Parameter a [yrs]  (the apparent initial age of the initial organic matter mixture ) as fct of sedimentation rate dependent formulation for a (Arndt et al., 2013)
                    load './data/Tmp_BW_hr_updated.mat'                             % BW temperature [°C]
                    load './data/sed_holo_hr_updated.mat'                           % Sedimentation rate [cm/yr] -- as in Bradley ea. 2020, after Burwicz ea. 2011
                    load './data/O2_BW_WOA2018_hr_updated.mat'                  	% BW  O2 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
                    load './data/NO3_BW_hr_updated.mat'                             % BW  NO3 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
                    load './data/FluxFe_total_Burw_BW_hr_updated.mat'               % Total FeOOH settling-flux [mumol/(cm2 yr)] -- just 16.67% of this is available for OM-degradation
                    load './data/water_depth_updated.mat'                           % Seafloor depth [m]
                    water_depth_updated = -water_depth_updated;
                    %          	load data/long.dat
                    load data/lat.dat
                    
                    % to calculate volume in grid-cell
                    dlat = 0.25;    % latitude difference in degrees for high-res
                    dlon = 0.25;    % longitude difference in degrees for high-res
                    
                case 2                                        
                    % load boundary conditions - all in 1 degree resolution
                    use_inversely_calculated_TOC = true;
                    use_Lee_toc_data = true;
                    if(use_Lee_toc_data)
                        if(use_inversely_calculated_TOC)
                            load('./data/01_plot_BC_April_21/BC_1degree/Lee_toc_lr_weighted_SWI_PK_filled_a1to10_por085.mat');  % inversely calculated TOC at SWI [wt%]     -- need to translate to mol/cm^3 solid phase, i.e. *1e-2/12*bsd.rho_sed
%                            load('./data/01_plot_BC_April_21/BC_1degree/Lee_toc_lr_weighted_SWI_a10to100_por085.mat');  % inversely calculated TOC at SWI [wt%]     -- need to translate to mol/cm^3 solid phase, i.e. *1e-2/12*bsd.rho_sed
                            toc_in = TOC_SWI_BC;
                        else
                            load('./data/01_plot_BC_April_21/BC_1degree/Lee_toc_lr_weighted.mat');       % Mean TOC in upper 5cm from Lee [wt%]     -- need to translate to mol/cm^3 solid phase, i.e. *1e-2/12*bsd.rho_sed
                            toc_in = Lee_toc_lr_weighted;
                        end
                        load('./data/01_plot_BC_April_21/BC_1degree/a_lr_aligned_PK_filled_in_a1to10.mat');                    %  Parameter a [yrs]  (the apparent initial age of the initial organic matter mixture ) as fct of sedimentation rate dependent formulation for a (Arndt et al., 2013)
                        a_hr_updated = Parameter_a_total_PKsedrate1to10;
%                        load('./data/01_plot_BC_April_21/BC_1degree/a_lr_aligned_10to100.mat');                               %  Parameter a [yrs]  (the apparent initial age of the initial organic matter mixture ) as fct of sedimentation rate dependent formulation for a (Arndt et al., 2013)
%                        load('./data/01_plot_BC_April_21/BC_1degree/a_lr_aligned_1to50.mat');                                 %  Parameter a [yrs]  (the apparent initial age of the initial organic matter mixture ) as fct of sedimentation rate dependent formulation for a (Arndt et al., 2013)
%                        a_hr_updated = a_lr_aligned;
                        load('./data/01_plot_BC_April_21/BC_1degree/sed_lr_weighted_aligned.mat');                           % Sedimentation rate [cm/yr] -- as in Bradley ea. 2020, after Burwicz ea. 2011
                        sed_holo_hr_updated = sed_lr_weighted_aligned;
                        load('./data/01_plot_BC_April_21/BC_1degree/water_depth_aligned.mat');                             % Seafloor depth [m]
                        water_depth_updated = -water_depth_aligned;
                        
                        load('./data/01_plot_BC_April_21/BC_1degree/lat_lr.mat');
                        lat = lat_lr;
                        
                        load('./data/01_plot_BC_April_21/BC_1degree/Tmp_BW_lr_aligned.mat');   % BW temperature [°C]
                        Tmp_BW_hr_updated = Tmp_BW_lr_aligned;
                        
                        load('./data/01_plot_BC_April_21/BC_1degree/O2_BW_WOA2018_lr_aligned.mat');                   	% BW  O2 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
                        O2_BW_WOA2018_hr_updated = O2_BW_WOA2018_lr_aligned;
                        load('./data/01_plot_BC_April_21/BC_1degree/NO3_BW_lr_aligned.mat');                             % BW  NO3 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
                        NO3_BW_hr_updated = NO3_BW_lr_aligned;
                        %                load './data/BC_1degrees/FluxFe_total_Burw_BW_01_updated.mat'               % Total FeOOH settling-flux [mumol/(cm2 yr)] -- just 16.67% of this is available for OM-degradation; used por=0.7 to calculate it
                        %                FluxFe_total_Burw_BW_hr_updated = FluxFe_total_Burw_BW_01_updated;
                        load('./data/01_plot_BC_April_21/BC_1degree/FluxFe_total_Burw_BW_aligned.mat'); % Total FeOOH settling-flux [mumol/(cm2 yr)] -- just 16.67% of this is available for OM-degradation; used por=0.85 to calculate it
                        FluxFe_total_Burw_BW_hr_updated = FluxFe_total_Burw_BW_01_por85_aligned;
                        %         	load data/BC_1degrees/long_WOA_01.mat
                        
                        
                        % Seiter_and_old_BC
                    else
                        load './data/BC_1degrees/toc_01.mat'                                       % TOC at SWI [wt%]     -- need to translate to mol/cm^3 solid phase, i.e. *1e-2/12*bsd.rho_sed
                        toc_in = toc_01;
                        load './data/BC_1degrees/a_01_updated.mat'                                  %  Parameter a [yrs]  (the apparent initial age of the initial organic matter mixture ) as fct of sedimentation rate dependent formulation for a (Arndt et al., 2013)
                        a_hr_updated = a_01_updated;
                        load './data/BC_1degrees/Tmp_BW_01_updated.mat'                             % BW temperature [°C]
                        Tmp_BW_hr_updated = Tmp_BW_01_updated;
                        load './data/BC_1degrees/sed_holo_01_updated.mat'                           % Sedimentation rate [cm/yr] -- as in Bradley ea. 2020, after Burwicz ea. 2011
                        sed_holo_hr_updated = sed_holo_01_updated;
                        load './data/BC_1degrees/O2_BW_WOA2018_01_updated.mat'                  	% BW  O2 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
                        O2_BW_WOA2018_hr_updated = O2_BW_WOA2018_01_updated;
                        load './data/BC_1degrees/NO3_BW_01_updated.mat'                             % BW  NO3 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
                        NO3_BW_hr_updated = NO3_BW_01_updated;
                        %                load './data/BC_1degrees/FluxFe_total_Burw_BW_01_updated.mat'               % Total FeOOH settling-flux [mumol/(cm2 yr)] -- just 16.67% of this is available for OM-degradation; used por=0.7 to calculate it
                        %                FluxFe_total_Burw_BW_hr_updated = FluxFe_total_Burw_BW_01_updated;
                        load './data/BC_1degrees/FluxFe_total_Burw_BW_01_por85_updated.mat'               % Total FeOOH settling-flux [mumol/(cm2 yr)] -- just 16.67% of this is available for OM-degradation; used por=0.85 to calculate it
                        FluxFe_total_Burw_BW_hr_updated = FluxFe_total_Burw_BW_01_por85_updated;
                        load './data/BC_1degrees/water_depth_updated.mat'                           % Seafloor depth [m]
                        water_depth_updated = -water_depth_updated;
                        %         	load data/BC_1degrees/long_WOA_01.mat
                        load data/BC_1degrees/lat_WOA_01.mat
                        lat = lat_WOA_01;
                    end
                    
                    % to calculate volume in grid-cell
                    dlat = 1.0;     % latitude difference in degrees for 1 degree resolution
                    dlon = 1.0;     % longitude difference in degrees for 1 degree resolution
                    
                case 3
                    % load boundary conditions - all in 2 degree resolution
                    use_inversely_calculated_TOC = true;
                    use_Lee_toc_data = true;
                    if(use_Lee_toc_data)
                        if(use_inversely_calculated_TOC)
                            load('./data/01_plot_BC_April_21/BC_2degree/Lee_TOC_SWI_BC_PK_and_SA_1to10_updated.mat');       % Inversely calculated TOC at SWI from mean toc in upper 5cm from Lee [wt%]     -- need to translate to mol/cm^3 solid phase, i.e. *1e-2/12*bsd.rho_sed
                            toc_in = Lee_TOC_SWI_BC_PK_and_SA_1to10_updated;
                        else
                            load('./data/01_plot_BC_April_21/BC_2degree/Lee_toc_lr_weighted.mat');       % Mean TOC in upper 5cm from Lee [wt%]     -- need to translate to mol/cm^3 solid phase, i.e. *1e-2/12*bsd.rho_sed
                            toc_in = Lee_toc_lr_weighted;
                        end
                        
%                         load('./data/01_plot_BC_April_21/BC_2degree/a_values_SA_1to100_updated.mat');                                   %  Parameter a [yrs]  (the apparent initial age of the initial organic matter mixture ) as fct of sedimentation rate dependent formulation for a (Arndt et al., 2013)
%                         a_hr_updated = a_values_SA_1to100_updated;
                        load('./data/01_plot_BC_April_21/BC_2degree/a_values_PK_and_SA_1to10_updated.mat');                                   %  Parameter a [yrs]  (the apparent initial age of the initial organic matter mixture ) as fct of sedimentation rate dependent formulation for a (Arndt et al., 2013)
                        a_hr_updated = a_values_PK_and_SA_1to10_updated;

                        load('./data/01_plot_BC_April_21/BC_2degree/sed_holo_02.mat');                           % Sedimentation rate [cm/yr] -- as in Bradley ea. 2020, after Burwicz ea. 2011
                        sed_holo_hr_updated = sed_holo_02;
                        load('./data/01_plot_BC_April_21/BC_2degree/water_depth_02_updated.mat');                             % Seafloor depth [m]
                        water_depth_updated = -water_depth_02_updated;
                        
                        load('./data/01_plot_BC_April_21/BC_2degree/lat_lr.mat');
                        lat = lat_lr;
                        
                        load('./data/01_plot_BC_April_21/BC_2degree/Tmp_BW_02_updated.mat');   % BW temperature [°C]
                        Tmp_BW_hr_updated = Tmp_BW_02_updated;
                        
                        load('./data/01_plot_BC_April_21/BC_2degree/O2_BW_WOA2018_02_updated.mat');                   	% BW  O2 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
                        O2_BW_WOA2018_hr_updated = O2_BW_WOA2018_02_updated;
                        load('./data/01_plot_BC_April_21/BC_2degree/NO3_BW_02_updated.mat');                             % BW  NO3 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
                        NO3_BW_hr_updated = NO3_BW_02_updated;
                        %                load './data/BC_1degrees/FluxFe_total_Burw_BW_01_updated.mat'               % Total FeOOH settling-flux [mumol/(cm2 yr)] -- just 16.67% of this is available for OM-degradation; used por=0.7 to calculate it
                        %                FluxFe_total_Burw_BW_hr_updated = FluxFe_total_Burw_BW_01_updated;
                        load('./data/01_plot_BC_April_21/BC_2degree/FluxFe_total_por085_Burw_BW_02_updated.mat'); % Total FeOOH settling-flux [mumol/(cm2 yr)] -- just 16.67% of this is available for OM-degradation; used por=0.85 to calculate it
                        FluxFe_total_Burw_BW_hr_updated = FluxFe_total_por085_Burw_BW_02_updated;
                        %         	load data/BC_1degrees/long_WOA_01.mat
                        
                        
                        % Seiter_and_old_BC

                    else
                    % load boundary conditions - all in 2 degree resolution
                    load './data/BC_2degrees/toc_02.mat'                                       % TOC at SWI [wt%]     -- need to translate to mol/cm^3 solid phase, i.e. *1e-2/12*bsd.rho_sed
                    toc_in = toc_02;
                    % a-values as fct of sedimentation rate:
                    load './data/BC_2degrees/a_02_updated.mat'                                  %  Parameter a [yrs]  (the apparent initial age of the initial organic matter mixture ) as fct of sedimentation rate dependent formulation for a (Arndt et al., 2013)
                    a_hr_updated = a_02_updated;
                    % Philip's a-values
%                     load './data/BC_2degrees/a_values_PK_updated.mat'                                  %  Parameter a [yrs]  (the apparent initial age of the initial organic matter mixture ) as fct of sedimentation rate dependent formulation for a (Arndt et al., 2013)
%                     a_hr_updated = a_values_PK_updated;
                    
                    load './data/BC_2degrees/Tmp_BW_02_updated.mat'                             % BW temperature [°C]
                    Tmp_BW_hr_updated = Tmp_BW_02_updated;
                    load './data/BC_2degrees/sed_holo_02.mat'                           % Sedimentation rate [cm/yr] -- as in Bradley ea. 2020, after Burwicz ea. 2011
                    sed_holo_hr_updated = sed_holo_02;
                    load './data/BC_2degrees/O2_BW_WOA2018_02_updated.mat'                  	% BW  O2 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
                    O2_BW_WOA2018_hr_updated = O2_BW_WOA2018_02_updated;
                    load './data/BC_2degrees/NO3_BW_02_updated.mat'                             % BW  NO3 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
                    NO3_BW_hr_updated = NO3_BW_02_updated;
                    %                load './data/BC_2degrees/FluxFe_total_Burw_BW_02_updated.mat'               % Total FeOOH settling-flux [mumol/(cm2 yr)] -- just 16.67% of this is available for OM-degradation; used por=0.7 to calculate it
                    %                FluxFe_total_Burw_BW_hr_updated = FluxFe_total_Burw_BW_02_updated;
                    load './data/BC_2degrees/FluxFe_total_Burw_BW_02_por85_updated.mat'               % Total FeOOH settling-flux [mumol/(cm2 yr)] -- just 16.67% of this is available for OM-degradation; used por=0.85 to calculate it
                    FluxFe_total_Burw_BW_hr_updated = FluxFe_total_Burw_BW_02_por85_updated;
                    load './data/BC_2degrees/water_depth_02_updated.mat'                           % Seafloor depth [m]
                    water_depth_updated = water_depth_02_updated;
                    %         	load data/BC_2degrees/long_WOA_02.mat
                    load data/BC_2degrees/lat_02.mat
                    %            long = long_WOA_02;
                    lat = lat_02;
            end
                    % to calculate volume in grid-cell
                    dlat = 2.0;     % latitude difference in degrees for 2 degree resolution
                    dlon = 2.0;     % longitude difference in degrees for 2 degree resolution
                    
            end
            
            if(OneGmodel)   % run as 1G with k-values as in Thullner ea. '09;
                matrix_k = zeros(size(a_hr_updated));
                matrix_k(water_depth_updated<150) = 0.221;
                matrix_k((water_depth_updated>=150) & (water_depth_updated<350)) = 0.208;
                matrix_k((water_depth_updated>=350) & (water_depth_updated<750)) = 0.174;
                matrix_k((water_depth_updated>=750) & (water_depth_updated<1500)) = 0.130;
                matrix_k((water_depth_updated>=1500) & (water_depth_updated<2750)) = 0.0718;
                matrix_k((water_depth_updated>=2750) & (water_depth_updated<4250)) = 0.0296;
                matrix_k(water_depth_updated>=4250) = 0.0122;
            else
                if(depth_dep_a)
                    a_hr_updated(water_depth_updated<50) = 1.0; % 3e-4;% 0.001; % shelf and slope environments after Dale and Thullner
                    a_hr_updated((water_depth_updated>=50) & (water_depth_updated<100)) = 10.0;%0.1; % shelf and slope environments after Dale and Thullner
                    %                a_hr_updated(water_depth_updated<100) = 3e-4;% 0.001; % shelf and slope environments after Dale and Thullner
                    a_hr_updated((water_depth_updated>=100) & (water_depth_updated<500)) = 50.0;%0.1; % shelf and slope environments after Dale and Thullner
                    a_hr_updated((water_depth_updated>=500) & (water_depth_updated<2000)) = 1.0;%0.1; % shelf and slope environments after Dale and Thullner
                    a_hr_updated((water_depth_updated>=2000) & (water_depth_updated< 3000)) = 10.0; % intermediate environments compare Bradley and Thullner
                    a_hr_updated((water_depth_updated>=3000) & (water_depth_updated< 4250)) = 25.0; % intermediate environments after Dale
                    a_hr_updated(water_depth_updated>=4250) = 100.0; % abyss
                else % use either Philip's a-values or sed-rate dependent but with some variations
                    
                    
                    %% this is not needed anymore when the a-values are already adjusted
                    rescale_a_values = false;                    
                    if(rescale_a_values)
                        a_hr_updated_in = a_hr_updated;
                        a_max = max(max(a_hr_updated));
                        a_min = min(min(a_hr_updated));
                        a_max_new = 10;
                        a_min_new = 1.0;
                        a_hr_updated = (a_hr_updated - a_min)./(a_max-a_min).*(a_max_new - a_min_new) + a_min_new;
                    end
                end
            end
            % counters for Problems in Fe-calculation
            debug.iter_fun6 = 0;                    % boundaries have same sign -> reduce gammaFe2
            debug.iter_Fefail = 1;                  % Problem macthing desired Fe3-influx (also zFe3 < zinf) start with 1 bc I initialized it with [0 0]
            debug.iter_zbioMatching_GMD = 0;        % counter for calculating the correct TOC profile using the GMD IntConst
           	debug.iter_zbioMatching_final = 0;      % counter for calculating the correct TOC profile after trying different methods
           	debug.iter_NEGATIVE_OXID = 0;           % counter for negative total OM oxidation

            [m,n]=size(toc_in);
            tic
            %% loop through lat - lon
            xstart = 1;  % 22
            ystart = 1;  % 65
            for x=xstart:m
                fprintf('-------------- New latitude ----------- \n');
                fprintf('x, %i \n',  x);
                
                % calculate volume
                % convert deg to cm CODAS package (by E.Firing,et al.) --
                % http://pordlabs.ucsd.edu/matlab/coord.htm --
                % lat - is the average latitude between the 2 fixes (lats)
                rlat = lat(x) * pi/180;
                m = 111132.09  - 566.05 * cos(2 * rlat)+ 1.2 * cos(4 * rlat);
                dy = dlat*m;        % latitude difference in meters
                p = 111415.13 * cos(rlat) - 94.55 * cos(3 * rlat);
                dx = dlon*p;        % longitude difference in meters
                
                for y=ystart:n
                    %                 fprintf('\n');
                    fprintf('x, y, %i %i \n',  x, y);
                    % now check for sed_holo or toc = NaN -- if yes set output to NaN
                    if ((isnan(sed_holo_hr_updated(x,y))) | (isnan(toc_in(x,y))) | (isnan(a_hr_updated(x,y))) | (isnan(Tmp_BW_hr_updated(x,y))) | (isnan(O2_BW_WOA2018_hr_updated(x,y))) | (isnan(water_depth_updated(x,y))))
                        dxdy(x,y)   = NaN;
                        
                        Cox_rate_out.Cox_aerobic(x,y) = NaN;
                        Cox_rate_out.Cox_denitr(x,y) = NaN;
                        Cox_rate_out.Cox_FeIII(x,y) = NaN;
                        Cox_rate_out.Cox_sulfate(x,y) = NaN;
                        Cox_rate_out.Cox_total(x,y) = NaN;
                        Cox_rate_out.Total_rate_sum(x,y) = NaN;
                        
                        Penetration_out.zox(x,y) = NaN;
                        Penetration_out.zno3(x,y) = NaN;
                        Penetration_out.zfeIII(x,y) = NaN;
                        Penetration_out.zso4(x,y) = NaN;
                        
                        SWI_fluxes.flxswiO2(x,y) = NaN;
                        if(swi.Nitrogen)
                            SWI_fluxes.flxswiNO3(x,y) = NaN;
                            SWI_fluxes.flxswiNH4(x,y) = NaN;
                        end
                        if(swi.Iron)
                            SWI_fluxes.flxswiFe2(x,y) = NaN;
                            Flux_Fe2_Dale_units(x,y) =  NaN;   	% in units of umol m-2 d-1 as in Dale ea. 2015, Fig. 3
                        end
                        SWI_fluxes.flxswiH2S(x,y) = NaN;
                        SWI_fluxes.flxswiSO4(x,y) = NaN;
                        
                        % save TOC flux (burial) at zinf in mol cm-2 yr-1
                        TOC_burial_flux(x,y) = NaN;
                        
                        
                        %                         % not necessary - results are essentially the
                        %                         same as whe using calReac above
                        %                         Cox_rate_out_flux.Cox_aerobic(x,y) = NaN;
                        %                         Cox_rate_out_flux.Cox_denitr(x,y) = NaN;
                        %                         Cox_rate_out_flux.Cox_FeIII(x,y) = NaN;
                        %                         Cox_rate_out_flux.Cox_sulfate(x,y) = NaN;
                        %                         Cox_rate_out_flux.Cox_Total(x,y) = NaN;
                        
                    else        % valid BC: run the model
                        

                        rho_sed_loc = 2.5;
                        swi.C0_nonbio = toc_in(x,y)*1e-2/12*rho_sed_loc;   % (wt%) -> (mol/cm^3 solid phase)
                        swi.p_a = a_hr_updated(x,y);
                        %                        swi.k = matrix_k(x,y);     % if use 1G model - but this leads to crazy OM degrdation
                        %swi.p_a = 15; %a_hr_updated(x,y);
                        swi.T = Tmp_BW_hr_updated(x,y);
                        swi.BC_sed_rate_flag = true;    % flag for existing BC for sed-rate
                        swi.BC_sed_rate = sed_holo_hr_updated(x,y);     % holocene sed-rate
                        swi.O20 = O2_BW_WOA2018_hr_updated(x,y)*10^(-9);
                        swi.NO30 = NO3_BW_hr_updated(x,y)*10^(-9);
                        Flux_FeIII0_in = Fe3_HR_frac*FluxFe_total_Burw_BW_hr_updated(x,y)*10^(-6);
                        if(swi.O20 > hypoxic_th)
                            swi.Flux_FeIII0 = (1-Fe3_in_SMI_frac_oxic)*Flux_FeIII0_in;  % just reactivre part in mol/(cm2 yr)
                            Flux_Fe2_from_SMI = Fe3_in_SMI_frac_oxic*Flux_FeIII0_in;  % in mol/(cm2 yr)

                        else
                            swi.Flux_FeIII0 = (1-Fe3_in_SMI_frac_hypoxic)*Flux_FeIII0_in;  % just reactivre part in mol/(cm2 yr)
                            Flux_Fe2_from_SMI = Fe3_in_SMI_frac_hypoxic*Flux_FeIII0_in;  % in mol/(cm2 yr)
                            Cox_rate_out.grid_hypoxic=Cox_rate_out.grid_hypoxic+1;
                            
                        end
                        swi.Flux_Fe2_from_SMI = Flux_Fe2_from_SMI;

                        %                  	Fe_influx_Dale = 0.1667*1110.0E-006;                     % FeIII influx from Dale but just 16.67% is available for DIR (See his Tab. 2, footnote d)
                        %                     swi.Flux_FeIII0 =  (Fe_influx_Dale)*365/100^2;           % Dale 1110 mol/(m^2 day)   -->  mol/(cm^2 yr):     *365/100^2
                        
                        swi.BC_wdepth_flag = true;
                        swi.BC_wdepth = water_depth_updated(x,y);
                        if(swi.BC_wdepth>7900)
                            swi.BC_wdepth %=7900;
                        end
                        
                        % RUN OMEN-SED
                        res=benthic_test.test_benthic(1,swi);
                        
                        % check also if we get a negative total oxidation rate -> in that case record as fail
                    	Cox_total_test = res.zTOC_RCM.calcReac(0.0, res.bsd.zinf, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2;
                        if(Cox_total_test<0.0)
                            res.zbio_Matching_fails = true;
                            fprintf('\n');
                            fprintf(' ---------  NEGATIVE TOTAL OXIDATION RATE AT x y %i %i \n',  x, y);
                            debug.iter_NEGATIVE_OXID = debug.iter_NEGATIVE_OXID+1;
                            NEGATIVE_OXID_xy(debug.iter_NEGATIVE_OXID, :) = [x y];                            
                        end                        
                        
                        % Did the calculation of the TOC profile fail -- or do we get a negative total OM oxidation?
                        if(res.zbio_Matching_fails)
                            fprintf('\n');
                            fprintf('zbio Matching with GMD IntConst fails at x y %i %i \n',  x, y);
                            debug.iter_zbioMatching_GMD = debug.iter_zbioMatching_GMD+1;
                            zbioMatching_xy_GMD(debug.iter_zbioMatching_GMD, :) = [x y];
                            
                            calculate_with_relaxed_zbioMatching = true;
                            if(calculate_with_relaxed_zbioMatching)     
                                % use SA Integration constants here they usually work 
                                % but are not 100% mass conserving
                            	swi.IntConst_GMD= false;         % true A1, A2 as in GMD; false: use Sandra's calculation
                                res=benthic_test.test_benthic(1,swi);
                            	swi.IntConst_GMD= true;         % true A1, A2 as in GMD; false: use Sandra's calculation         
                                
                                % check again for a negative total oxidation rate -> in that case record as fail
                                Cox_total_test = res.zTOC_RCM.calcReac(0.0, res.bsd.zinf, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2;
                                if(Cox_total_test<0.0)
                                    res.zbio_Matching_fails = true;
                                    fprintf('\n');
                                    fprintf('------   STILL NEGATIVE TOTAL OXIDATION RATE AT x y %i %i \n',  x, y);
                                    debug.iter_NEGATIVE_OXID = debug.iter_NEGATIVE_OXID+1;
                                    NEGATIVE_OXID_xy(debug.iter_NEGATIVE_OXID, :) = [x y];                            
                                end 

                            end
                        end
                        
                        % check again for failed matching (or negative
                        % total Cox, in case it is negative with the new
                        % Integration Constants
                        if(res.zbio_Matching_fails) % TOC matching failed - save locations as NaN
                            fprintf('\n');
                            fprintf('zbio Matching fails at x y %i %i \n',  x, y);
                            debug.iter_zbioMatching_final = debug.iter_zbioMatching_final+1;
                            zbioMatching_xy_final(debug.iter_zbioMatching_final, :) = [x y];
                            
                            % Set all result-variables to NaN
                            dxdy(x,y)   = NaN;
                            
                            Cox_rate_out.Cox_aerobic(x,y) = NaN;
                            Cox_rate_out.Cox_denitr(x,y) = NaN;
                            Cox_rate_out.Cox_FeIII(x,y) = NaN;
                            Cox_rate_out.Cox_sulfate(x,y) = NaN;
                            Cox_rate_out.Cox_total(x,y) = NaN;
                            Cox_rate_out.Total_rate_sum(x,y) = NaN;
                            
                            Penetration_out.zox(x,y) = NaN;
                            Penetration_out.zno3(x,y) = NaN;
                            Penetration_out.zfeIII(x,y) = NaN;
                            Penetration_out.zso4(x,y) = NaN;
                            
                            SWI_fluxes.flxswiO2(x,y) = NaN;
                            if(swi.Nitrogen)
                                SWI_fluxes.flxswiNO3(x,y) = NaN;
                                SWI_fluxes.flxswiNH4(x,y) = NaN;
                            end
                            if(swi.Iron)
                                SWI_fluxes.flxswiFe2(x,y) = NaN;
                                Flux_Fe2_Dale_units(x,y) =  NaN;   	% in units of umol m-2 d-1 as in Dale ea. 2015, Fig. 3
                            end
                            SWI_fluxes.flxswiH2S(x,y) = NaN;
                            SWI_fluxes.flxswiSO4(x,y) = NaN;
                            
                            % save TOC flux (burial) at zinf in mol cm-2 yr-1
                            TOC_burial_flux(x,y) = NaN;
                            
                        else
                            debug.iter_fun6 = debug.iter_fun6 + res.iter_fun6;
                            debug.iter_Fefail = debug.iter_Fefail + res.iter_Fefail;
                            % save location for fauiled Fe-calculation
                            if(res.iter_Fefail > 0)
                                Fe_fail_xy(debug.iter_Fefail, :) = [x y];
                            end
                            
                            dxdy(x,y)   = dx*dy;    % m^2 per grid-cell
                            
                            
                            swi.C0i = res.swi.C0i;
                            % oxidation rates in mmol m-2 yr-1
                            Cox_rate_out.Cox_aerobic(x,y) = res.zTOC_RCM.calcReac(0.0, res.zox, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2;
                            if(swi.Nitrogen)
                                Cox_rate_out.Cox_denitr(x,y) = res.zTOC_RCM.calcReac(res.zox, res.zno3, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2;
                            end
                            if(res.swi.Iron)
                                Cox_rate_out.Cox_FeIII(x,y) = res.zTOC_RCM.calcReac(res.zno3, res.zfeIII,1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2;
                            end
                            Cox_rate_out.Cox_sulfate(x,y) = res.zTOC_RCM.calcReac(res.zfeIII,res.zso4, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2;
                            Cox_rate_out.Cox_total(x,y) = res.zTOC_RCM.calcReac(0.0, res.bsd.zinf, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2;
                            Cox_rate_out.Total_rate_sum(x,y) = Cox_rate_out.Cox_aerobic(x,y) + Cox_rate_out.Cox_denitr(x,y) + Cox_rate_out.Cox_FeIII(x,y) + Cox_rate_out.Cox_sulfate(x,y);
                            
                            Penetration_out.zox(x,y) = res.zox;
                            Penetration_out.zno3(x,y) = res.zno3;
                            Penetration_out.zfeIII(x,y) = res.zfeIII;
                            Penetration_out.zso4(x,y) = res.zso4;
                            
                            % SWI-fluxes in mol cm-2 yr-1
                            SWI_fluxes.flxswiO2(x,y) = res.flxswiO2;
                            if(swi.Nitrogen)
                                SWI_fluxes.flxswiNO3(x,y) = res.flxswiNO3;
                                SWI_fluxes.flxswiNH4(x,y) = res.flxswiNH4;
                            end
                            if(swi.Iron)
                                SWI_fluxes.flxswiFe2(x,y) = res.flxswiFe2;
                                Flux_Fe2_Dale_units(x,y) = res.flxswiFe2 *10^6 *100^2 /365;   	% in units of umol m-2 d-1 as in Dale ea. 2015, Fig. 3
                            end
                            SWI_fluxes.flxswiH2S(x,y) = res.flxswiH2S;
                            SWI_fluxes.flxswiSO4(x,y) = res.flxswiSO4;
                            
                            % compare results with Dale
                            [F_TOC_swi, F_TOC1_swi] = res.zTOC_RCM.calcCflx(0, res.bsd, res.swi, res);
                            [F_TOC_inf, F_TOC1_inf] = res.zTOC_RCM.calcCflx(res.bsd.zinf, res.bsd, res.swi, res);
                            % save TOC flux (burial) at zinf in mol cm-2 yr-1
                            TOC_burial_flux(x,y) = F_TOC_inf;
                            
                            Cox_total_diff = (F_TOC_swi-F_TOC_inf)*1000 *100^2/365;      % in units of mmol m-2 d-1 as in Dale ea. 2015, Fig. 2a
                            
                            
                            %                         % calculate mean OM concentration in upper x cm
                            %                         depth = 5.0;
                            %                         [C_swi, C1_swi] = res.zTOC_RCM.calcC( 0.0, res.bsd, res.swi, res);
                            %                         [C_swi2, C1_swi2] = res.zTOC_RCM.calcC( 0.1, res.bsd, res.swi, res);
                            %                         [C_5, C1_5] = res.zTOC_RCM.calcC( depth, res.bsd, res.swi, res);
                            %                         [C_95, C1_95] = res.zTOC_RCM.calcC( 9.5, res.bsd, res.swi, res);
                            %                         [C_zbio, C1_zbio] = res.zTOC_RCM.calcC( 10.0, res.bsd, res.swi, res);
                            %                         [C_zbio_pl, C1_zbio_pl] = res.zTOC_RCM.calcC( 10.0, res.bsd, res.swi, res);
                            %                         OM_0cm=C_swi* 100*12/res.bsd.rho_sed;
                            %                         OM_01cm=C_swi2* 100*12/res.bsd.rho_sed;
                            %                         OM_5cm=C_5* 100*12/res.bsd.rho_sed;
                            %                         OM_95cm=C_95* 100*12/res.bsd.rho_sed;
                            %                         OM_zbio=C_zbio* 100*12/res.bsd.rho_sed;
                            %                         OM_zbio_pl=C_zbio_pl* 100*12/res.bsd.rho_sed;
                            %                         Mean_OM = 1/depth * 100*12/res.bsd.rho_sed*res.zTOC_RCM.calcOM(0.0, depth, 1, res.bsd, res.swi, res);
                            %
                            %                         str_date = [datestr(date,7), datestr(date,5), datestr(date,11)];
                            %                         benthic_test.plot_column(res, false, swi, str_date)
                            
                        end     % Did TOC matching at zbio fail?
                    end
                end
            end
            toc
            %% plot output
            %            benthic_test.plot_Fe_SA(swi.POC_flux, O20, Cox_rate_out, Penetration_out, Flux_Fe2_Dale_units, swi.nG)
            
            str_date = [datestr(date,7), datestr(date,5), datestr(date,11)];
            
            save(['./output/' exp_name '_Cox_rate_out_' str_date '.mat'] , 'Cox_rate_out')
            save(['./output/' exp_name '_Flux_Fe2_Dale_units_' str_date '.mat'] , 'Flux_Fe2_Dale_units')
            save(['./output/' exp_name '_Penetration_out_' str_date '.mat'] , 'Penetration_out')
            save(['./output/' exp_name '_Fe_fail_xy_' str_date '.mat'] , 'Fe_fail_xy')
            %            save(['./output/Cox_rate_out_flux' exp_name '_' str_date '.mat'] , 'Cox_rate_out_flux')
            save(['./output/' exp_name '_dxdy_' str_date '.mat'] , 'dxdy')
            save(['./output/' exp_name '_SWI_fluxes_' str_date '.mat'] , 'SWI_fluxes')
            save(['./output/' exp_name '_TOC_burial_flux_' str_date '.mat'] , 'TOC_burial_flux')
            
            debug.Fe_fail_xy = Fe_fail_xy;
            debug.zbioMatching_xy_GMD = zbioMatching_xy_GMD;
            debug.zbioMatching_xy_final = zbioMatching_xy_final;
            debug.NEGATIVE_OXID_xy = NEGATIVE_OXID_xy;
        end
        
        
        function [Cox_rate_out, Penetration_out, Flux_Fe2_Dale_units] = run_OMEN_RCM_testDale_POC(ExpName)
            % run OMEN-SED with default SWI conditions as in default_swi()

            swi = benthic_test.default_swi();
            swi.TwoG_OM_model = false;
            swi.Test_Dale = true;
            swi.Test_Dale_14G = true;
            
            Fe3_HR_frac = 0.1667;                                           % 16.67% of total FeOOH-influx is HR and thus available for DIR (See his Tab. 2, footnote d)
            Flux_FeIII0_in = Fe3_HR_frac*1110.0E-006*365/100^2;           % Dale 1110 mol/(m^2 day)   -->  mol/(cm^2 yr):     *365/100^2
            
            hypoxic_th = 20.0e-9;           % weakly hypoxic threshold for hypoxia after Dale et al. '15
            anoxic_th = 3.0e-9;           % weakly hypoxic threshold for hypoxia after Dale et al. '15
            Fe3_in_SMI_frac_oxic = 0.6;     % Fe3_HR used for Sulphide-Mediated Iron reduction under oxic BW (i.e. not available for DIR)          
            Fe3_in_SMI_frac_hypoxic = 0.28;  % Fe3_HR used for Sulphide-Mediated Iron reduction under hypoxic BW (i.e. not available for DIR)
            Fe3_in_SMI_frac_anoxic = 0.23;  % Fe3_HR used for Sulphide-Mediated Iron reduction under hypoxic BW (i.e. not available for DIR)
            
            if(swi.Test_Dale)
                swi.BC_wdepth_flag = true;
                swi.BC_wdepth = 350;    % for Dale experiment
                swi.BC_sed_rate_flag = true;
                swi.BC_sed_rate = 60/1000;  % To compare with Dale et al. (2015)
                if(swi.Test_Dale_14G)
                    swi.nG = 14;
                end
                swi.p_a = 10.0;
                swi.p_nu = 0.125;
                
            end
            
            swi.sed_column_depth = 30;
            swi.POC_flux = [0.5 1 2 4 6 8 10 12 14 16];
            swi.flux = true;
            O20 = [1e-9 2e-9 5e-9 10e-9 15e-9 25e-9 50e-9 100e-9 200e-9];
            
            for i=1:length(swi.POC_flux)
                i
                for j=1:length(O20)
                    fprintf('\n');
                    fprintf('i, j %i %i \n',  i, j);
                    
                    swi.POCi = i;
                    swi.O20 = O20(j);
                    if(swi.O20 > hypoxic_th)
                        swi.Flux_FeIII0 = (1-Fe3_in_SMI_frac_oxic)*Flux_FeIII0_in;  % just reactivre part in mol/(cm2 yr)
                        Flux_Fe2_from_SMI = Fe3_in_SMI_frac_oxic*Flux_FeIII0_in;  % in mol/(cm2 yr)
                    elseif(swi.O20 > anoxic_th)
                        	swi.Flux_FeIII0 = (1-Fe3_in_SMI_frac_hypoxic)*Flux_FeIII0_in;  % just reactivre part in mol/(cm2 yr)
                            Flux_Fe2_from_SMI = Fe3_in_SMI_frac_hypoxic*Flux_FeIII0_in;  % in mol/(cm2 yr)
                    else
                        swi.Flux_FeIII0 = (1-Fe3_in_SMI_frac_anoxic)*Flux_FeIII0_in;  % just reactivre part in mol/(cm2 yr)
                    	Flux_Fe2_from_SMI = Fe3_in_SMI_frac_anoxic*Flux_FeIII0_in;  % in mol/(cm2 yr)
                    end
                    swi.Flux_Fe2_from_SMI = Flux_Fe2_from_SMI;
                    
                    res=benthic_test.test_benthic(1,swi);
                    
                    Flux_Fe2_from_SMI_Daleunits = Flux_Fe2_from_SMI*10^6*100^2/365;     % from mol/(cm2 yr) to umol m-2 d-1
                    % set date-time or string going into plot function
                    
                    swi.C0i = res.swi.C0i;
                    Cox_rate_out.Cox_aerobic(i,j) = res.zTOC_RCM.calcReac(0.0, res.zox, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2/365;
                    if(swi.Nitrogen)
                        Cox_rate_out.Cox_denitr(i,j) = res.zTOC_RCM.calcReac(res.zox, res.zno3, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2/365;
                    end
                    if(res.swi.Iron)
                        if(res.zno3 == res.bsd.zinf)
                            Cox_rate_out.Cox_FeIII(i,j) = 0.0;
                        else
                            Cox_rate_out.Cox_FeIII(i,j) = res.zTOC_RCM.calcReac(res.zno3, res.zfeIII,1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2/365;
                        end
                    end
                    Cox_rate_out.Cox_sulfate(i,j) = res.zTOC_RCM.calcReac(res.zfeIII,res.zso4, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2/365;
                    Cox_rate_out.Cox_total(i,j) = res.zTOC_RCM.calcReac(0.0, res.bsd.zinf, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2/365;
                    Cox_rate_out.Total_rate_sum(i,j) = Cox_rate_out.Cox_aerobic(i,j) + Cox_rate_out.Cox_denitr(i,j) + Cox_rate_out.Cox_FeIII(i,j) + Cox_rate_out.Cox_sulfate(i,j);
                    
                    Cox_frac_Fe_red(i,j)=Cox_rate_out.Cox_FeIII(i,j)/Cox_rate_out.Cox_total(i,j)*100;
                    Cox_frac_Fe_red_sum(i,j)=Cox_rate_out.Cox_FeIII(i,j)/Cox_rate_out.Total_rate_sum(i,j)*100;
                    
                    Penetration_out.zox(i,j) = res.zox;
                    Penetration_out.zno3(i,j) = res.zno3;
                    Penetration_out.zfeIII(i,j) = res.zfeIII;
                    Penetration_out.zso4(i,j) = res.zso4;
                    % compare results with Dale
                    [F_TOC_swi, F_TOC1_swi] = res.zTOC_RCM.calcCflx(0, res.bsd, res.swi, res);
                    [F_TOC_inf, F_TOC1_inf] = res.zTOC_RCM.calcCflx(res.bsd.zinf, res.bsd, res.swi, res);
                    
                    Cox_total_diff = (F_TOC_swi-F_TOC_inf)*1000 *100^2/365;      % in units of mmol m-2 d-1 as in Dale ea. 2015, Fig. 2a
                    Flux_Fe2_Dale_units(i,j) = res.flxswiFe2 *10^6 *100^2 /365;                   % in units of umol m-2 d-1 as in Dale ea. 2015, Fig. 3
%                    Flux_Fe2_from_SMI_Daleunits_clean = (1-res.bsd.gammaFe_pp)*(1 - res.bsd.gammaFe2)*Flux_Fe2_from_SMI_Daleunits;
%                    Flux_Fe2_Dale_units(i,j) = res.flxswiFe2 *10^6 *100^2 /365 + Flux_Fe2_from_SMI_Daleunits_clean;
                    
                    
                    %% old output:
                    % calculate degradation pathways
                    %                 Cox_rate_out(1,2*i-1) = Cox_rate.Cox_aerobic;
                    %                 Cox_rate_out(2,2*i-1) = Cox_rate.Cox_denitr;
                    %                 Cox_rate_out(3,2*i-1) = Cox_rate.Cox_FeIII;
                    %                 Cox_rate_out(4,2*i-1) = Cox_rate.Cox_sulfate;
                    %                 Cox_rate_out(5,2*i-1) = Cox_rate.Cox_total;
                    %                 % FRaction of total:
                    %                 Cox_rate_out(1,2*i) = Cox_rate.Cox_aerobic/Cox_rate.Cox_total*100;
                    %                 Cox_rate_out(2,2*i) = Cox_rate.Cox_denitr/Cox_rate.Cox_total*100;
                    %                 Cox_rate_out(3,2*i) = Cox_rate.Cox_FeIII/Cox_rate.Cox_total*100;
                    %                 Cox_rate_out(4,2*i) = Cox_rate.Cox_sulfate/Cox_rate.Cox_total*100;
                    %                 Cox_rate_out(5,2*i) = Cox_rate.Cox_total/Cox_rate.Cox_total*100;
                    %
                    %
                    %                 Penetration_out(1,2*i-1) = res.zox;
                    %                 Penetration_out(2,2*i-1) = res.zno3;
                    %                 Penetration_out(3,2*i-1) = res.zfeIII;
                    %                 Penetration_out(4,2*i-1) = res.zso4;
                    %                 Penetration_out(1,2*i) = nan;
                    %                 Penetration_out(2,2*i) = nan;
                    %                 Penetration_out(3,2*i) = nan;
                    %                 Penetration_out(4,2*i) = nan;
                end
            end
            
            %% plot output
            benthic_test.plot_Fe_SA(swi.POC_flux, O20, Cox_rate_out, Penetration_out, Flux_Fe2_Dale_units, swi.nG, swi.sed_column_depth, swi.p_a, ExpName)
            
            
        end
        
        function res = test_benthic( ncl, swi )
            %            loc_BW_O2_anoxia = 5.0e-9;       	% set to 5.0 nanomol/cm^3
            
            % counters for problems in calculation
            res.iter_fun6 = 0;      % boundaries have same sign -> reduce gammaFe2
            res.iter_Fefail = 0;    % Problem macthing desired Fe3-influx (also zFe3 < zinf)
            
            if nargin < 1
                ncl = 1;
            end
            
            if(swi.BC_wdepth_flag)
                if(swi.BC_wdepth<=200) % change porosity for shelfs, i.e. <200m
                    por_shelf = 0.85;
                	res.bsd = benthic_main(ncl, swi.BC_wdepth, por_shelf);   % specify 350m water-depth
                else
                    res.bsd = benthic_main(ncl, swi.BC_wdepth);   % specify 350m water-depth
                end
                res.bsd.zinf = swi.sed_column_depth;
                

            else
                res.bsd = benthic_main(ncl);
            end
            res.bsd.usescalarcode = ncl==1;
            
            if(swi.BC_sed_rate_flag)
                res.bsd.w = swi.BC_sed_rate;               
                %                 if(res.bsd.w>0.1)
                %                 res.bsd.w = 0.0005;    % To compare with Dale et al. (2015)
                %                 end
            end
            
            if(swi.Test_Dale)
               % change bioturbation coefficient (and zbio) as in Dale et al. (2015)
               a = 20;  % μM
               b = 12;  % μM
               O2BW = swi.O20*10^9;      %  μM
               f = 0.5+0.5.*erf((O2BW-a)/b);               
               res.bsd.Dbio = 23.0 * f;
               res.bsd.zbio = 7.0 * f;
               
            end
            
            if nargin < 2 || isempty(swi)
                swi = benthic_test.default_swi();
            end
            
            if ncl > 1  % set up O2 gradient for testing
                O20 = swi.O20;
                for i = 1:ncl
                    swi.O20(i) = 10*(i-1)/(ncl-1)*O20;
                end
            end
            
            res.swi = swi;
            
            % initialize then calculate
            if(swi.TwoG_OM_model)
                res.zTOC = benthic_zTOC(res.bsd);
            else
                res.zTOC_RCM = benthic_zTOC_RCM(res.bsd);
            end
            
            %            tic;
            if(swi.TwoG_OM_model)
                res = res.zTOC.calc(res.bsd,res.swi, res);
                O2_demand_flux = -(res.swi.Fnonbio1+res.swi.Fnonbio2)*res.bsd.OC/((1-res.bsd.por)./res.bsd.por);
                res.zbio_Matching_fails = false;
%                 % plot 2G TOC profile?
%                 benthic_test.plot_TOC_2G(res, false, res.swi, '0706_2G_OMEN')

            else
                % Adding into on RCM for MultiG approach
                if(res.swi.Test_Dale)
                    [res.zTOC_RCM.k, res.swi.C0i, res.swi.Fnonbioi] = benthic_test.RCM_Dale_2015(res.bsd, res.swi);
                else
                    [res.zTOC_RCM.k, res.swi.C0i, res.swi.Fnonbioi] = benthic_test.RCM(res.bsd, res.swi);
                end
                res = res.zTOC_RCM.calc(res.bsd,res.swi, res);
                O2_demand_flux = -(sum(res.swi.Fnonbioi))*res.bsd.OC/((1-res.bsd.por)./res.bsd.por);
            end
            
            % % Check TOC matching at zbio
            % If is is not correct: give NaNs back
            if(res.zbio_Matching_fails)
                % Don't do any calculations here. The TOC profile is wrong!
                
            else
                
                if(swi.TwoG_OM_model)
                    %    Cox_total = res.zTOC.calcReac(0.0, res.bsd.zinf, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2/365; % calculate in mmol/(m2 d) for Seb's parameteriztion
                    
                else
                    Cox_total = res.zTOC_RCM.calcReac(0.0, res.bsd.zinf, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2/365; % calculate in mmol/(m2 d) for Seb's parameteriztion
                end
                res.zO2 = benthic_zO2(res.bsd, res.swi);
                res.zNO3 = benthic_zNO3(res.bsd, res.swi);
                res.zFeIII = benthic_zFeIII(res.bsd, res.swi);
                res.zSO4 = benthic_zSO4(res.bsd, res.swi);
                res.zNH4 = benthic_zNH4(res.bsd, res.swi);
                res.zFe2 = benthic_zFe2(res.bsd, res.swi, Cox_total);
                res.zH2S = benthic_zH2S(res.bsd, res.swi);
                if(swi.calc_P_DIC_ALK)
                    res.zPO4_M = benthic_zPO4_M(res.bsd, res.swi);
                    res.zDIC = benthic_zDIC(res.bsd, res.swi);
                    res.zALK = benthic_zALK(res.bsd, res.swi);
                end
                % set gamma for Fe as fct of Cox and O20:
                BW_O2 = swi.O20 * 1000 *10^6;               % converst from mol/cm3 to mol/kg = mol/L (*1000), then express as mu (*10^6)
                % calculate gamma for fraction of Fe2 re-oxidation at zox, using Seb vd Velde's fit
                a=0.89;
                b=-3.83;
                c=1.45;
                x = log(Cox_total/BW_O2);
                gammaFe2_max = a*exp(-((x-b)^2)/(2*c^2));
                % gammaFe2_max = 0.0;
                % save fractions up to this value in vector in case
                % Fe-calculation fails and I need to iteratively decrease gammaFe2
                gammaFe2_v = [0:0.1:gammaFe2_max gammaFe2_max];
                gammaFe2_v = fliplr(gammaFe2_v);
                % calculate gamma for Fe2 burial (mainly as pyrite), using Seb vd Velde's fit
                d=0.154 * 6;                                % because FFeOOh_t = 6 * FFeOOH_hr
                res.bsd.gammaFe_pp = (1-d*tanh(Cox_total/BW_O2));       % 1 - parameterization for DFe-efflux
                %            res.bsd.gammaFe_pp = 0.5;
                
                
                % %     change zbio and gammaH2S depending on BW oxygenation:
                % %             if(res.swi.O20<=loc_BW_O2_anoxia)   % anoxic sediments?
                % %                 bsd.zbio = 0.01;        % decrease bioturbation depth
                % %                 bsd.gammaH2S = 0.95;    % fraction of H2S that is oxidised in anoxic sediments
                % %             end
                
                % loop over gammaFe2:
                % in case Fe-cycle is not mass-conserving gammaFe2 has to be smaller
                % this just happens for a few cases when zno3 < zbio << zfe3
                % (but I am not sure why) - calculation leads to zfe3 < zinf
                % but not all Fe3flux_in is used
                for(iter=1:length(gammaFe2_v))           % iteratively chose smaller gammaFe2
                    %                iter
                    res.bsd.gammaFe2 = gammaFe2_v(iter);
                    if(res.swi.O20<=0.0)
                        res.zox=0.0;
                        res.flxzox = 0.0;
                        res.conczox = 0.0;
                        res.flxswiO2=0.0;
                        res.zxf=0.0;
                    else
                        res = res.zO2.calc(res.bsd, res.swi, res);
                    end
                    if(swi.Nitrogen)
                        if(res.swi.NO30<=0.0)
                            res.zno3=res.zox;
                            res.flxzno3 = 0.0;
                            res.conczno3 = 0.0;
                            res.flxswiNO3 = 0.0;
                        else
                            res = res.zNO3.calc(res.bsd, res.swi, res);
                        end
                    else
                        res.zno3=res.zox;
                        %                res.zso4=res.zox;   % for test-case with just TOC & O2
                    end
                    
                    if(swi.Iron)
                        % First solve iteratively for correct SWI [Fe3+]
                        Flux_FeIII0_in = 2;
                        fun=@(Fe3_conc_in) -res.zFeIII.calc_input(res.bsd, res.swi, res, Fe3_conc_in) - res.swi.Flux_FeIII0;   % NOTE: negative calc_input to transform the result to positive because the given flux at the SWI is positive here!
                        
                        % For starting point calculate SWI-concentration of Fe3+ from input flux
                        if(res.zno3 > res.bsd.zbio )
                            % SWI-concentration of Fe3+ is not affected by Fe-reduction
                            % -> use simple conversion - this should be close
                            % enough to correct SWI-concentration
                            res.swi.FeIII0= res.swi.Flux_FeIII0/((1-res.bsd.por)*res.bsd.w);     % calculate concentration [mol/cm^3] from flux [mol/(cm2 yr)] according non-bioturbated flux!!!
                            %                    test0_flx = res.zFeIII.calc_input(res.bsd, res.swi, res, res.swi.FeIII0);
                            % try zero flux at zinf and see if we have any FeIII left, also
                            % calculate [FeIII] at zinf for advective loss
                            [flxzfeIII, conczinf, flxswi,rtmp] = res.zFeIII.calcbc(res.bsd.zinf, res.bsd, res.swi, res, 2);
                            if(conczinf >= 0)    % zfeIII >= zinf --> use first guess for fe3-flux and assume the remainder is buried
                                fun0 = 0.0;     % use the first guess for res.swi.FeIII0
                                % FOR MASS BALANCE WITH cGENIE MAKE SURE THE
                                % REMAINING INFLUX OF Fe3+ IS BURIED IN SEDIMENTS!!!
                            else    % zfeIII < zinf --> iteratively solve for correct Fe3-concentration
                                fun0 = fun(res.swi.FeIII0);     % initial test with calculated SWI-concentration
                            end
                        else
                            % SWI-concentration of Fe3+ is affected by Fe-reduction
                            % -> take biodiffusion and Fe-reduction into account
                            res.swi.FeIII0 = res.zFeIII.calcConcentration(res.bsd, res.swi, res, res.swi.Flux_FeIII0);
                            if(res.swi.FeIII0 < 0)
                                % Either Case zno3 < zFe3 < zbio OR zno3 < zbio << zFe3 (= zinf) and calcConcentration() fails
                                % use simple conversion and iterate from there
                                res.swi.FeIII0= res.swi.Flux_FeIII0/((1-res.bsd.por)*res.bsd.w);     % calculate concentration [mol/cm^3] from flux [mol/(cm2 yr)] according non-bioturbated flux!!!
                                [flxzfeIII, conczinf, flxswi,rtmp] = res.zFeIII.calcbc(res.bsd.zinf, res.bsd, res.swi, res, 2);
                                % check if calculated SWI-flux for calculated concentration and 'infinite' concentration
                                % is smaller than desired SWI-flux (i.e. check for res.swi.FeIII0 & 1e+6);
                                %                        fun_FeIII00 = fun(1e-10);
                                fun_FeIII0 = fun(res.swi.FeIII0);
                                %                         fun_FeIIIinf = fun(1e+6); % DH deleted this check 02.11.2020, problem converging if zfeIII = zinf???
                                % If yes: zfeIII = inf bc SWI-flux of fe3 won't increase anymore if SWI-conc. is increased
                                % -> use first guess for fe3-flux
                                %                        if(conczinf >= 0)    % zfeIII >= zinf --> use first guess for fe3-flux and assume the remainder is buried
                                %                        if((fun_FeIII0 < 0) & (fun_FeIIIinf < 0))    % DH deleted the 2nd check 02.11.2020, problem converging if zfeIII = zinf???
                                if((fun_FeIII0 < 0))% & (fun_FeIIIinf < 0))    % zbio < zfeIII < zinf --> use first guess for fe3-flux and assume the remainder is buried
                                    fun0 = 0.0;     % use the first guess for res.swi.FeIII0
                                else    % zfeIII < zinf --> iteratively solve for correct Fe3-concentration
                                    %                         test0_flx = res.zFeIII.calc_input(res.bsd, res.swi, res, res.swi.FeIII0);
                                    fun0 = fun(res.swi.FeIII0);     % initial test with calculated SWI-concentration
                                end
                                %                          fun0 = fun(res.swi.FeIII0);     % initial test with calculated SWI-concentration
                            else
                                % try zero flux at zinf and see if we have any FeIII left
                                [flxzfeIII, conczinf, flxswi,rtmp] = res.zFeIII.calcbc(res.bsd.zinf, res.bsd, res.swi, res, 2);
                                if(conczinf >= 0)    % zfeIII >= zinf --> use first guess for fe3-flux and assume the remainder is buried
                                    fun0 = 0.0;     % use the first guess for res.swi.FeIII0
                                else    % zfeIII < zinf --> iteratively solve for correct Fe3-concentration
                                    % check if calculated SWI-flux for calculated concentration and 'infinite' concentration
                                    % is smaller than desired SWI-flux (i.e. check for res.swi.FeIII0 & 1e+6);
                                    %                        fun_FeIII00 = fun(1e-10);
                                    fun_FeIII0 = fun(res.swi.FeIII0);
                                    %                        fun_FeIIIinf = fun(1e+6);  % DH deleted this check 02.11.2020, problem converging if zfeIII = zinf???
                                    % If yes: zfeIII = inf bc SWI-flux of fe3 won't increase anymore if SWI-conc. is increased
                                    % -> use first guess for fe3-flux
                                    %                        if((fun_FeIII0 < 0) & (fun_FeIIIinf < 0))    % DH deleted the 2nd check 02.11.2020, problem converging if zfeIII = zinf???
                                    if((fun_FeIII0 < 0))    % zbio < zfeIII < zinf --> no Fe2-reoxidation
                                        fun0 = 0;
                                    else    % zfeIII < zinf --> iteratively solve for correct Fe3-concentration
                                        %                         test0_flx = res.zFeIII.calc_input(res.bsd, res.swi, res, res.swi.FeIII0);
                                        fun0 = fun_FeIII0; %fun(res.swi.FeIII0);     % initial test with calculated SWI-concentration
                                    end
                                    %                            fun0 = fun(res.swi.FeIII0);     % initial test with calculated SWI-concentration
                                end
                            end
                        end
                        
                        if(abs(fun0) > res.bsd.tol_Fe3)  % initial guess above not good enough
                            
                            if(fun0 > 0)    % Fe3+ concentration too high
                                test0_flx = fun(1e-10);    % res.zFeIII.calc_input(res.bsd, res.swi, res, 1e-10);
                                test1_flx = fun(res.swi.FeIII0);         % res.zFeIII.calc_input(res.bsd, res.swi, res, res.swi.FeIII0);
                                res.swi.FeIII0=fzero(fun,[1e-10, res.swi.FeIII0],res.bsd.fzerooptions_Fe3);    % NOTE: Find better way for lower boundary
                            else    % Fe3+ concentration too low
                                %                                                   fun2=fun(2.0*res.swi.FeIII0);
                                fun6 = fun(1);
                                %                                                  test3_flx = res.zFeIII.calc_input(res.bsd, res.swi, res, 1e+6);
                                %                                                  test4_flx = res.zFeIII.calc_input(res.bsd, res.swi, res, 1e+3);
                                if(fun6<0)    % check if fun for upper and lower boundary both below 0 --> use lower gammaFe2 value!
                                    % do nothing just use first guess for res.swi.FeIII0
                                    % if not good enough reduce gammaFe2 in next iteration
                                    res.iter_fun6 =  res.iter_fun6+1;
                                else
                                    res.swi.FeIII0=fzero(fun,[res.swi.FeIII0, 1.0],res.bsd.fzerooptions_Fe3);    % now just use 20% higher FeIII0 input
                                    % was until 02.12.2020  res.swi.FeIII0=fzero(fun,[res.swi.FeIII0, 1e+6],res.bsd.fzerooptions_Fe3);    % NOTE: Find better way for upper boundary - not sure if this works anyway (at least not for gammaFe2 >= 0.99)
                                end
                            end
                        end
                        res = res.zFeIII.calc(res.bsd, res.swi, res);
                        % check IF FeIII-flux_in is calculated or IF zfe3 = zinf
                        % otherwise simulate with gammaFe2 = 0
                        if((abs(-res.flxswiFeIII-res.swi.Flux_FeIII0) < res.bsd.tol_Fe3) || (res.zfeIII == res.bsd.zinf))
                            % all good
                            break    % leave for-loop
                            %                    res.flxswiFeIII
                            %                    res.zfeIII
                        else
                            warning('SOMETHING WRONG WITH Fe-FLUXES!');
                            % reduce re-oxidation of Fe2 , in next for loop
                            
                            if(iter==length(gammaFe2_v))     % does not work for gammaFe2_v = 0
                                warning('FeIII calculation fails for gammaFe2_v = 0')
                                res.iter_Fefail =  res.iter_Fefail+1;
                            end
                            
                        end
                    else
                        res.zfeIII = res.zno3;
                        break    % leave for-loop
                    end
                end
                
                res = res.zSO4.calc(res.bsd, res.swi, res);
                if(swi.Nitrogen)
                    res = res.zNH4.calc(res.bsd, res.swi, res);
                end
                if(swi.Iron)
                    
                    res = res.zFe2.calc(res.bsd, res.swi, res);
                    % correct Fe2 return flux by Fe2 lost to pyrite
                    % precipitation and SMI
                    res.flxswiFe2 = res.flxswiFe2 + swi.Flux_Fe2_from_SMI*(1-res.bsd.gammaFe_pp)*(1 - res.bsd.gammaFe2);
                end
                res = res.zH2S.calc(res.bsd, res.swi, res);
                if(swi.calc_P_DIC_ALK)
                    res = res.zPO4_M.calc(res.bsd, res.swi, res);
                    res = res.zDIC.calc(res.bsd, res.swi, res);
                    res = res.zALK.calc(res.bsd, res.swi, res);
                end
                %            toc;
                
                %%%%% WRITE OUTPUT:
                if(swi.write_output)
                    answ = res
                    if(swi.TwoG_OM_model)
                        % use fluxes to calculate burial efficiency and total Cox:
                        [F_TOC_swi, F_TOC1_swi] = res.zTOC.calcCflx(0, res.bsd, res.swi, res);
                        [F_TOC_inf, F_TOC1_inf] = res.zTOC.calcCflx(res.bsd.zinf, res.bsd, res.swi, res);
                        
                        [F_TOC_zox, F_TOC1_zox] = res.zTOC.calcCflx(res.zox, res.bsd, res.swi, res);
                        [F_TOC_zno3, F_TOC1_zno3] = res.zTOC.calcCflx(res.zno3, res.bsd, res.swi, res);
                        [F_TOC_zfe3, F_TOC1_zfe3] = res.zTOC.calcCflx(res.zfeIII, res.bsd, res.swi, res);
                        [F_TOC_zso4, F_TOC1_zso4] = res.zTOC.calcCflx(res.zso4, res.bsd, res.swi, res);
                        
                        fprintf('Fluxes at zinf %g \n',  F_TOC_inf);
                        fprintf('Fluxes at swi %g \n',  F_TOC_swi);
                        
                        fprintf('Burial efficiency of POC (in %%) %g \n',  F_TOC_inf/F_TOC_swi *100);
                        fprintf('\n');
                        fprintf('Total Cox (mol cm-2 yr-1) %e \n',  F_TOC_swi-F_TOC_inf);
                        fprintf('Aerobic Cox %e \n',  F_TOC_swi-F_TOC_zox);
                        fprintf('Denitrification Cox %e \n',  F_TOC_zox - F_TOC_zno3 );
                        fprintf('Fe-reduction Cox %e \n',  F_TOC_zno3 - F_TOC_zfe3);
                        fprintf('Sulfate reduction Cox %e \n',  F_TOC_zfe3-F_TOC_zso4);
                        
                    else
                        % use fluxes to calculate burial efficiency and total Cox:
                        [F_TOC_swi, F_TOC1_swi] = res.zTOC_RCM.calcCflx(0, res.bsd, res.swi, res);
                        [F_TOC_inf, F_TOC1_inf] = res.zTOC_RCM.calcCflx(res.bsd.zinf, res.bsd, res.swi, res);
                        
                        [F_TOC_zox, F_TOC1_zox] = res.zTOC_RCM.calcCflx(res.zox, res.bsd, res.swi, res);
                        [F_TOC_zno3, F_TOC1_zno3] = res.zTOC_RCM.calcCflx(res.zno3, res.bsd, res.swi, res);
                        [F_TOC_zfe3, F_TOC1_zfe3] = res.zTOC_RCM.calcCflx(res.zfeIII, res.bsd, res.swi, res);
                        [F_TOC_zso4, F_TOC1_zso4] = res.zTOC_RCM.calcCflx(res.zso4, res.bsd, res.swi, res);
                        
                        fprintf('Fluxes at zinf %g \n',  F_TOC_inf);
                        fprintf('Fluxes at swi %g \n',  F_TOC_swi);
                        fprintf('Burial efficiency of POC (in %%) %g \n',  F_TOC_inf/F_TOC_swi *100);
                        fprintf('\n');
                        fprintf('Total Cox (mol cm-2 yr-1) %e \n',  F_TOC_swi-F_TOC_inf);
                        fprintf('Aerobic Cox %e \n',  F_TOC_swi-F_TOC_zox);
                        fprintf('Denitrification Cox %e \n',  F_TOC_zox - F_TOC_zno3 );
                        fprintf('Fe-reduction Cox %e \n',  F_TOC_zno3 - F_TOC_zfe3);
                        fprintf('Sulfate reduction Cox %e \n',  F_TOC_zfe3-F_TOC_zso4);
                        
                        %             %%% WRITE EXACT FLUX
                        %             FO2_exact=res.zO2.calcFO2_exact(res.zox,res.bsd, res.swi, res);
                        %             fprintf('exact flxswiO2 flux (mol cm^{-2} yr^{-1}) %g \n',  FO2_exact);
                        
                    end
                end
                
            end     % Check TOC matching at zbio
        end
        
        
        
        function res = calculate_meanTOC( ncl, swi )
            %% set OMEN and then only calculate mean TOC in upper 5cm
            
            if nargin < 1
                ncl = 1;
            end
            
            if(swi.BC_wdepth_flag)
                res.bsd = benthic_main(ncl, swi.BC_wdepth);   % specify 350m water-depth
                res.bsd.zinf = swi.sed_column_depth;
            else
                res.bsd = benthic_main(ncl);
            end
            res.bsd.usescalarcode = ncl==1;
            
            if(swi.BC_sed_rate_flag)
                res.bsd.w = swi.BC_sed_rate;
            end
            
            if nargin < 2 || isempty(swi)
                swi = benthic_test.default_swi();
            end
            
            res.swi = swi;
            
            % initialize then calculate
            res.zTOC_RCM = benthic_zTOC_RCM(res.bsd);
            
            [res.zTOC_RCM.k, res.swi.C0i, res.swi.Fnonbioi] = benthic_test.RCM(res.bsd, res.swi);
            res = res.zTOC_RCM.calc(res.bsd,res.swi, res);
            
            % calculate mean OM concentration in upper x cm
            x = 5;
            [C_x, C1_11] = res.zTOC_RCM.calcC( x, res.bsd, res.swi, res);
            OM_10=C_x* 100*12/res.bsd.rho_sed;     % OM wt% at x cm
            res.Mean_OM = 1/x * 100*12/res.bsd.rho_sed*res.zTOC_RCM.calcOM(0.0, x, 1, res.bsd, swi, res);   % mean OM et% in upper 5cm
            
            
        end
        
        function plot_Fe_SA(POC_flux, O20, Cox_rate_out, Penetration_out, Flux_Fe2_Dale_units, nG, column_depth, a_param, ExpName)
            %% plot various results for the Fe analysis
            
            str_date = datestr(now,'ddmmyy');
            vs_Cox = true;
            plot_effluxes = true;
            plot_frac_pathways = false;
            
            if(vs_Cox)
                plot_name = '_vs_Cox';
                txt_xLabel = '[Cox] (mmol m^{-2} d^{-1})';
            else                
                plot_name = '_vs_POCflx';
                txt_xLabel = 'POCflx (mmol m^{-2} d^{-1})';
                % make a matrix from POC-influx for plotting
                POC_tr= POC_flux';
                POCflx_matrix = repmat(POC_tr,1,length(O20));

            end
            
            if(plot_effluxes)
            %% plot FFe efflux
            x_axis = [0 14];
            y_axis = [0 200];
            
            % plot FFe2 vs BW O2
            fig1 = figure('Renderer', 'painters', 'Position', [10 10 600 300]);
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            box on
            color = parula(length(POC_flux));
            hold on
            for i=1:length(POC_flux)
                plot(O20, Flux_Fe2_Dale_units(i,:),'x-','Color',color(i,:))
            end
            ylim(y_axis)
            ylabel('JDFe (\mumol m^{-2} d^{-1})');
            xlabel('[O_2]_{BW} (\muM)');
            
            print(fig1, '-depsc2', ['Flux_JFFe2_vs_O2_OMEN_' num2str(nG) 'G_' num2str(a_param) 'a_' num2str(column_depth) 'cm_' ExpName , '_' str_date '.eps']);
            
            % plot FFe2 vs Cox
            fig2 = figure('Renderer', 'painters', 'Position', [10 10 600 300]);
            set(0,'defaultLineLineWidth', 2)
            box on
            set(0,'DefaultAxesFontSize',12)
            color = parula(length(O20));
            hold on
            for j=1:length(O20)
                if(vs_Cox)
                    plot(Cox_rate_out.Cox_total(:,j), Flux_Fe2_Dale_units(:,j),'x-','Color',color(j,:))
                else
                    plot(POCflx_matrix(:,j), Flux_Fe2_Dale_units(:,j),'x-','Color',color(j,:))
                end
            end
%             xlim(x_axis)
            ylim(y_axis)
            ylabel('JDFe (\mumol m^{-2} d^{-1})');
            xlabel(txt_xLabel);
            
            print(fig2, '-depsc2', ['Flux_JFFe2' plot_name '_' num2str(nG) 'G_' num2str(a_param) 'a_' num2str(column_depth) 'cm_' ExpName , '_' str_date '.eps']);
            end
            
            if(plot_frac_pathways)
            %% Plot fraction of metabolic pathways
            y_axis = [0 100];
            % plot fraction aerobic-reduction vs Cox
            % plot FFe2 vs Cox
            fig3 = figure('Renderer', 'painters', 'Position', [10 10 600 300]);
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            box on
            color = parula(length(O20));
            hold on
            for j=1:length(O20)
                if(vs_Cox)
                    plot(Cox_rate_out.Cox_total(:,j), Cox_rate_out.Cox_aerobic(:,j)./Cox_rate_out.Cox_total(:,j)*100,'x-','Color',color(j,:))
                else
                    plot(POCflx_matrix(:,j), Cox_rate_out.Cox_aerobic(:,j)./Cox_rate_out.Cox_total(:,j)*100,'x-','Color',color(j,:))
                end
            end
            ylabel('Fract. of aerobic-reduction (%)');
            xlabel(txt_xLabel);
            %             xlim(x_axis)
            ylim(y_axis)
            print(fig3, '-depsc2', ['Frac_aerobic_red' plot_name '_' num2str(nG) 'G_' num2str(a_param) 'a_' num2str(column_depth) 'cm_' ExpName , '_' str_date  '.eps']);
            
            % plot fraction Fe-reduction vs Cox
            % plot FFe2 vs Cox
            fig4 = figure('Renderer', 'painters', 'Position', [10 10 600 300]);
            box on
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            color = parula(length(O20));
            hold on
            for j=1:length(O20)
                if(vs_Cox)
                    plot(Cox_rate_out.Cox_total(:,j), Cox_rate_out.Cox_denitr(:,j)./Cox_rate_out.Cox_total(:,j)*100,'x-','Color',color(j,:))
                else
                    plot(POCflx_matrix(:,j), Cox_rate_out.Cox_denitr(:,j)./Cox_rate_out.Cox_total(:,j)*100,'x-','Color',color(j,:))
                end
            end
            ylabel('Fract. of denitrification (%)');
            xlabel(txt_xLabel);
            %             xlim(x_axis)
            ylim(y_axis)
            print(fig4, '-depsc2', ['Frac_denitrif' plot_name '_' num2str(nG) 'G_'  num2str(a_param) 'a_' num2str(column_depth) 'cm_' ExpName , '_' str_date '.eps']);
            
            % plot fraction Fe-reduction vs Cox
            fig5 = figure('Renderer', 'painters', 'Position', [10 10 600 300]);
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            box on
            color = parula(length(O20));
            hold on
            for j=1:length(O20)
                if(vs_Cox)
                    plot(Cox_rate_out.Cox_total(:,j), Cox_rate_out.Cox_FeIII(:,j)./Cox_rate_out.Cox_total(:,j)*100,'x-','Color',color(j,:))
                else
                    plot(POCflx_matrix(:,j), Cox_rate_out.Cox_FeIII(:,j)./Cox_rate_out.Cox_total(:,j)*100,'x-','Color',color(j,:))
                end
            end
            ylabel('Fract. of Fe-reduction (%)');
            xlabel(txt_xLabel);
            %             xlim(x_axis)
           ylim([0 20])
            print(fig5, '-depsc2', ['Frac_Fe2_red' plot_name '_' num2str(nG) 'G_'  num2str(a_param) 'a_' num2str(column_depth) 'cm_' ExpName , '_' str_date '_zoom.eps']);
            
            % plot fraction Fe-reduction vs Cox
            % plot FFe2 vs Cox
            fig51 = figure('Renderer', 'painters', 'Position', [10 10 600 300]);
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            box on
            color = parula(length(O20));
            hold on
            for j=1:length(O20)
                if(vs_Cox)
                    plot(Cox_rate_out.Cox_total(:,j), Cox_rate_out.Cox_FeIII(:,j)./Cox_rate_out.Cox_total(:,j)*100,'x-','Color',color(j,:))
                else
                    plot(POCflx_matrix(:,j), Cox_rate_out.Cox_FeIII(:,j)./Cox_rate_out.Cox_total(:,j)*100,'x-','Color',color(j,:))
                end
            end
            ylabel('Fract. of Fe-reduction (%)');
            xlabel(txt_xLabel);
            ylim(y_axis)
            print(fig51, '-depsc2', ['Frac_Fe2_red' plot_name '_' num2str(nG) 'G_'  num2str(a_param) 'a_' num2str(column_depth) 'cm_' ExpName , '_' str_date  '.eps']);
            
            % plot fraction SO4-reduction vs Cox
            fig6 = figure('Renderer', 'painters', 'Position', [10 10 600 300]);
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            box on
            color = parula(length(O20));
            hold on
            for j=1:length(O20)
                if(vs_Cox)
                    plot(Cox_rate_out.Cox_total(:,j), Cox_rate_out.Cox_sulfate(:,j)./Cox_rate_out.Cox_total(:,j)*100,'x-','Color',color(j,:))
                else
                    plot(POCflx_matrix(:,j), Cox_rate_out.Cox_sulfate(:,j)./Cox_rate_out.Cox_total(:,j)*100,'x-','Color',color(j,:))
                end
            end
            ylabel('Fract. of SO_4-reduction (%)');
            xlabel(txt_xLabel);
            ylim(y_axis)
            %           xlim(x_axis)
            print(fig6, '-depsc2', ['Frac_SO4_red' plot_name '_' num2str(nG) 'G_'  num2str(a_param) 'a_' num2str(column_depth) 'cm_' ExpName , '_' str_date '.eps']);
            end
        end
        
        
        function plot_TOC_2G(res, debug, swi, str_date)
            % plot single sediment column vs depth
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            bsd = res.bsd;
            zgrid = 0:0.1:bsd.zinf;
            
            
            % CONCENTRATIONS WITHOUT PO4
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
                        for i=1:length(zgrid)
                            [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                            [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                        end
                        
                        color = parula(2);

                        %%% TOC wt %
                        fig_toc = figure
                        subplot(1,2,1)
                        plot(100*C1*12/bsd.rho_sed, -zgrid,'Color',color(2,:))
                        hold on
                        plot(100*C2*12/bsd.rho_sed, -zgrid,'Color',color(1,:))
                        plot(100*C*12/bsd.rho_sed, -zgrid, 'k')
                        t=xlim;         % to draw penetration depths the correct lengths
                        plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                        
                        %            plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                        hold off
                        yticks([-100:20:0])
                        xlabel ('TOC (wt%)')
                        ylabel('Depth (cm)')
                        %            title('Total TOC (wt%)')
                      	print(fig_toc, '-depsc2', [str_date '_TOC_PROFILES_2G.eps']);

                    end
        
        function plot_column(res, debug, swi, str_date)
            % plot single sediment column vs depth
            print_Fig1 = true;
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            bsd = res.bsd;
            zgrid = 0:0.1:bsd.zinf;
            
            
            % CONCENTRATIONS WITHOUT PO4
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            % Nitrogen included?
            if(print_Fig1)
                if(swi.Nitrogen)
                    %figure;
                    fig1 = figure('Renderer', 'painters', 'Position', [10 10 600 900]);
                    % TOC
                    if(swi.TwoG_OM_model)
                        subplot(3,2,1)
                        for i=1:length(zgrid)
                            [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                            [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                        end
                        
                        %%% TOC wt %
                        plot(100*C1*12/bsd.rho_sed, -zgrid, 'b')
                        hold on
                        plot(100*C2*12/bsd.rho_sed, -zgrid, 'g')
                        plot(100*C*12/bsd.rho_sed, -zgrid, 'k')
                        t=xlim;         % to draw penetration depths the correct lengths
                        plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                        plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                        plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                        plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                        
                        %            plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                        hold off
                        %                ylim([-50 0.0])
                        xlabel ('TOC (wt%)')
                        ylabel('Depth (cm)')
                        %            title('Total TOC (wt%)')
                    else
                        subplot(3,2,1)
                        for i=1:length(zgrid)
                            [C(i), C1(i,:)] = res.zTOC_RCM.calcC( zgrid(i), bsd, res.swi, res);
                            [Cflx(i), C1flx(i,:)] = res.zTOC_RCM.calcCflx( zgrid(i), bsd, res.swi, res);
                        end
                        %                 color = parula(swi.nG);
                        %                 for G = 1:swi.nG
                        %                     plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
                        %                     hold on
                        %                 end
                        % Plot sum (TOC)
                        plot(100*C*12/bsd.rho_sed, -zgrid, 'b')
                        hold on
                        t=xlim;         % to draw penetration depths the correct lengths
                        plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                        plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                        plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                        plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                        
                        hold off
                        xlabel ('TOC (wt%)')
                        ylabel('Depth (cm)')
                        %            title('Total TOC (wt%)')
                        
                    end
                    
                    %%% O2
                    if(res.zox>0.0)
                        for i=1:length(zgrid)
                            [O2(i), flxO2(i), flxO2D(i), flxO2adv(i)] = res.zO2.calcO2(zgrid(i), bsd, res.swi, res);
                        end
                    else
                        O2(i) = 0.0;
                    end
                    subplot(3,2,3)
                    plot(O2, -zgrid, 'b')
                    hold on
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    %                ylim([-20 0.0])
                    xlabel ('O_2 (mol/cm^3)')
                    ylabel('Depth (cm)')
                    %            title ('O2 (mol/cm^3)')
                    
                    %%% NO3
                    for i=1:length(zgrid)
                        [NO3(i), flxNO3(i)] = res.zNO3.calcNO3(zgrid(i), bsd, res.swi, res);
                    end
                    subplot(3,2,5)
                    plot(NO3, -zgrid, 'b')
                    hold on
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    %                ylim([-20 0.0])
                    xlabel ('NO_3 (mol/cm^3)')
                    ylabel('Depth (cm)')
                    %            title ('NO3 (mol/cm^3)')
                    
                    %%% NH4
                    for i=1:length(zgrid)
                        [NH4(i), flxNH4(i)] = res.zNH4.calcNH4(zgrid(i), bsd, res.swi, res);
                    end
                    subplot(3,2,4)
                    plot(NH4, -zgrid, 'b')
                    hold on
                    %                xlim([0 1e-7])
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    hold off
                    xlabel ('NH_4 (mol/cm^3)')
                    %                ylabel('Depth (cm)')
                    %            title ('NH4 (mol/cm^3)')
                    
                    %%% SO4
                    subplot(3,2,2)
                    for i=1:length(zgrid)
                        [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), bsd, res.swi, res);
                    end
                    plot(SO4, -zgrid, 'b')
                    hold on
                    %                xlim([2.7e-5 swi.SO40])
                    %                xlim([2.7e-5 swi.SO40])
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    hold off
                    xlabel ('SO_4 (mol/cm^3)')
                    %                ylabel('Depth (cm)')
                    %            title ('SO4 (mol/cm^3)')
                    
                    %%% H2S
                    subplot(3,2,6)
                    for i=1:length(zgrid)
                        [H2S(i), flxH2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                    end
                    plot(H2S, -zgrid, 'b')
                    hold on
                    %                xlim([0 4e-7])
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    xlabel ('H_2S (mol/cm^3)')
                    %               ylabel('Depth (cm)')
                    %            title ('H2S (mol/cm^3)')
                    
                    % save Figure
                    print(fig1, '-depsc2', [str_date '_ALL_PROFILES.eps']);
                    
                else % no Nitrogen
                    
                    figure;
                    % TOC
                    if(swi.TwoG_OM_model)
                        subplot(2,2,1)
                        for i=1:length(zgrid)
                            [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                            [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                        end
                        %%% TOC wt %
                        plot(100*C1*12/bsd.rho_sed, -zgrid, 'b')
                        hold on
                        plot(100*C2*12/bsd.rho_sed, -zgrid, 'g')
                        plot(100*C*12/bsd.rho_sed, -zgrid, 'k')
                        t=xlim;         % to draw penetration depths the correct lengths
                        plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                        plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                        plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                        plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                        
                        %            plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                        hold off
                        xlabel ('TOC (wt%)')
                        ylabel('Depth (cm)')
                        %            title('Total TOC (wt%)')
                    else
                        subplot(2,2,1)
                        for i=1:length(zgrid)
                            [C(i), C1(i,:)] = res.zTOC_RCM.calcC( zgrid(i), bsd, res.swi, res);
                            [Cflx(i), C1flx(i,:)] = res.zTOC_RCM.calcCflx( zgrid(i), bsd, res.swi, res);
                        end
                        color = parula(swi.nG);
                        for G = 1:swi.nG
                            plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
                            hold on
                        end
                        % Plot sum (TOC)
                        plot(100*C*12/bsd.rho_sed, -zgrid, ':k')
                        hold on
                        t=xlim;         % to draw penetration depths the correct lengths
                        plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                        plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                        plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                        plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                        
                        hold off
                        xlabel ('TOC (wt%)')
                        ylabel('Depth (cm)')
                        %            title('Total TOC (wt%)')
                        
                    end
                    
                    %%% O2
                    for i=1:length(zgrid)
                        [O2(i), flxO2(i), flxO2D(i), flxO2adv(i)] = res.zO2.calcO2(zgrid(i), bsd, res.swi, res);
                    end
                    subplot(2,2,3)
                    plot(O2, -zgrid, 'b')
                    hold on
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    xlabel ('O_2 (mol/cm^3)')
                    ylabel('Depth (cm)')
                    %            title ('O2 (mol/cm^3)')
                    
                    %%% SO4
                    subplot(2,2,2)
                    for i=1:length(zgrid)
                        [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), bsd, res.swi, res);
                    end
                    plot(SO4, -zgrid, 'b')
                    hold on
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    hold off
                    %            xlim([2.7e-5 swi.SO40])
                    xlabel ('SO_4 (mol/cm^3)')
                    ylabel('Depth (cm)')
                    %            title ('SO4 (mol/cm^3)')
                    
                    %%% H2S
                    subplot(2,2,4)
                    for i=1:length(zgrid)
                        [H2S(i), flxH2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                    end
                    plot(H2S, -zgrid, 'b')
                    hold on
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    xlabel ('H_2S (mol/cm^3)')
                    ylabel('Depth (cm)')
                    %            title ('H2S (mol/cm^3)')
                    
                    
                    print('-dpsc2', [str_date '_ALL_PROFILES.eps']);
                end
            end
            
            % Also plot the 2nd figure for PO4, DIC and ALK?
            if(swi.plot_PO4_DIC_ALK)
                % figure
                fig2 = figure('Renderer', 'painters', 'Position', [10 10 600 900]);
                
    % ACTUALLY DO NOT PLOT Fe-profiles as they are not correct! The
    % penetration depth and SWI-fluxes are okay though
%                 if(swi.Iron)
%                     %%% FeIII
%                     subplot(3,2,1)
%                     for i=1:length(zgrid)
%                         [FeIII(i), flxFeIII(i)] = res.zFeIII.calcFeIII(zgrid(i), bsd, res.swi, res);
%                     end
%                     plot(FeIII, -zgrid, 'b')
%                     hold on
%                     %                xlim([2.7e-5 swi.FeIII0])
%                     %                xlim([2.7e-5 swi.FeIII0])
%                     t=xlim;         % to draw penetration depths the correct lengths
%                     plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
%                     plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
%                     plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
%                     plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
%                     plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
%                     hold off
%                     %                ylim_min = -round(res.zfeIII,-1)-10;
%                     %                ylim([ylim_min 0.0])
%                     xlabel ('FeIII (mol/cm^3)')
%                     ylabel('Depth (cm)')
%                     %            title ('SO4 (mol/cm^3)')
%                     
%                     %%% Fe2
%                     subplot(3,2,2)
%                     for i=1:length(zgrid)
%                         [Fe2(i), flxFe2(i)] = res.zFe2.calcFe2(zgrid(i), bsd, res.swi, res);
%                     end
%                     plot(Fe2, -zgrid, 'b')
%                     hold on
%                     t=xlim;         % to draw penetration depths the correct lengths
%                     plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
%                     plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
%                     plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
%                     plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
%                     plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
%                     xlabel ('Fe2 (mol/cm^3)')
%                     %                ylabel('Depth (cm)')
%                     %            title ('Fe2 (mol/cm^3)')
%                     
%                 end
                %%% PO4
                subplot(3,2,3)
                for i=1:length(zgrid)
                    [PO4(i), flxPO4(i), M(i), flxM(i), e_M(i), f_M(i), p_M(i), q_M(i), g_M(i), dedz_M(i), dfdz_M(i), dpdz_M(i), dqdz_M(i), dgdz_M(i)] = res.zPO4_M.calcPO4_M(zgrid(i), bsd, res.swi, res);
                end
                plot(PO4, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                %          axis([0 1.5*10^(-9) -100 0])
                xlabel ('PO_4 (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('PO_4 (mol/cm^3)')
                
                %%% Fe-bound P (M)
                subplot(3,2,4)
                %for i=1:length(zgrid)
                %    [PO4(i), flxPO4(i), M(i), flxM(i)] = res.zPO4_M.calcPO4_M(zgrid(i), bsd, res.swi, res);
                %end
                plot(M, -zgrid, 'b')
                hold on
                %            plot([0,max(M)], [-bsd.zbio,-bsd.zbio], 'k--')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('Fe-bound P (mol/cm^3)')
                %            ylabel('Depth (cm)')
                %            title ('Fe-bound P (mol/cm^3)')
                
                %%% DIC
                subplot(3,2,5)
                for i=1:length(zgrid)
                    [DIC(i), flxDIC(i)] = res.zDIC.calcDIC(zgrid(i), bsd, res.swi, res);
                end
                plot(DIC, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('DIC (mol/cm^3)')
                ylabel('Depth (cm)')
                
                %%% ALK
                subplot(3,2,6)
                for i=1:length(zgrid)
                    [ALK(i), flxALK(i)] = res.zALK.calcALK(zgrid(i), bsd, res.swi, res);
                end
                plot(ALK, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('ALK (mol/cm^3)')
                ylabel('Depth (cm)')
                
                
                print(fig2, '-depsc2', [str_date '_PO4_PROFILES.eps']);
            end
            
            
            %%% some plots for debugging!
            if debug
                figure
                
                subplot(3,2,1)
                hold on
                plot(e_M, -zgrid, 'b')
                title ('Fe-P - ODE solution')
                
                subplot(3,2,2)
                hold on
                plot(f_M, -zgrid, 'b')
                
                subplot(3,2,3)
                hold on
                plot(p_M, -zgrid, 'b')
                
                subplot(3,2,4)
                hold on
                plot(q_M, -zgrid, 'b')
                
                subplot(3,2,5)
                hold on
                plot(g_M, -zgrid, 'b')
                
                figure
                hold on
                subplot(3,2,1)
                plot(dedz_M, -zgrid, 'b')
                title ('Fe-P - ODE derivations')
                subplot(3,2,2)
                plot(dfdz_M, -zgrid, 'b')
                subplot(3,2,3)
                plot(dpdz_M, -zgrid, 'b')
                subplot(3,2,4)
                plot(dqdz_M, -zgrid, 'b')
                subplot(3,2,5)
                plot(dgdz_M, -zgrid, 'b')
                
                
                %%%%%%%%%%%%%%%%%%%%%
                
                %         H2S
                
                %%%%%%%%%%%%%%%%%%%%%
                figure
                subplot(3,3,2)
                for i=1:length(zgrid)
                    [H2S(i), flxH2S(i), e_H2S(i), dedz_H2S(i), f_H2S(i), dfdz_H2S(i), g_H2S(i), dgdz_H2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                    % no func. called calcH2S_debug [H2S(i), flxH2S(i), e_H2S(i), dedz_H2S(i), f_H2S(i), dfdz_H2S(i), g_H2S(i), dgdz_H2S(i)] = res.zH2S.calcH2S_debug(zgrid(i), bsd, res.swi, res);
                    
                end
                plot(H2S, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('H_2S (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('H2S (mol/cm^3)')
                
                subplot(3,3,4)
                hold on
                plot(e_H2S, -zgrid, 'b')
                title ('H_2S - ODE solution')
                
                subplot(3,3,5)
                hold on
                plot(f_H2S, -zgrid, 'b')
                
                subplot(3,3,6)
                hold on
                plot(g_H2S, -zgrid, 'b')
                
                subplot(3,3,7)
                hold on
                plot(dedz_H2S, -zgrid, 'b')
                title ('H_2S - ODE derivations')
                
                subplot(3,3,8)
                hold on
                plot(dfdz_H2S, -zgrid, 'b')
                
                subplot(3,3,9)
                hold on
                plot(dgdz_H2S, -zgrid, 'b')
                
                
                
                % CONCENTRATION + Vertical Trransport
                figure;
                % TOC
                if(swi.TwoG_OM_model)
                    subplot(3,4,1)
                    for i=1:length(zgrid)
                        [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                        [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                    end
                    % TOC wt %
                    plot(100*C1*12/bsd.rho_sed, -zgrid, 'b')
                    hold on
                    plot(100*C2*12/bsd.rho_sed, -zgrid, 'g')
                    plot(100*C*12/bsd.rho_sed, -zgrid, 'k')
                    plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                    hold off
                    xlabel ('TOC (wt%)')
                    ylabel('Depth (cm)')
                    title('Total TOC (wt%)')
                    % TOC vertical transport flux
                    subplot(3,4,2);
                    plot(C1flx, -zgrid, 'b')
                    hold on
                    plot(C2flx, -zgrid, 'g')
                    plot(Cflx, -zgrid, 'k')
                    xlabel ('TOC trspt (mol cm^{-2}yr^{-1})')
                    ylabel('Depth (cm)')
                    title('TOC vert transport')
                else
                    subplot(3,4,1)
                    for i=1:length(zgrid)
                        [C(i), C1(i,:)] = res.zTOC_RCM.calcC( zgrid(i), bsd, res.swi, res);
                        [Cflx(i), C1flx(i,:)] = res.zTOC_RCM.calcCflx( zgrid(i), bsd, res.swi, res);
                    end
                    color = parula(swi.nG);
                    for G = 1:swi.nG
                        plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
                        hold on
                    end
                    % Plot sum (TOC)
                    plot(100*C*12/bsd.rho_sed, -zgrid, ':k')
                    hold on
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    hold off
                    xlabel ('TOC (wt%)')
                    ylabel('Depth (cm)')
                    % TOC vertical transport flux
                    subplot(3,4,2);
                    plot(C1flx, -zgrid, 'b')
                    hold on
                    plot(Cflx, -zgrid, 'k')
                    xlabel ('TOC trspt (mol cm^{-2}yr^{-1})')
                    ylabel('Depth (cm)')
                    title('TOC vert transport')
                    
                end
                
                % O2
                for i=1:length(zgrid)
                    [O2(i), flxO2(i), flxO2D(i), flxO2adv(i)] = res.zO2.calcO2(zgrid(i), bsd, res.swi, res);
                end
                subplot(3,4,3)
                plot(O2, -zgrid, 'b')
                hold on
                plot([0,res.swi.O20], [-bsd.zbio,-bsd.zbio], 'k--')
                xlabel ('O2 (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('O2 (mol/cm^3)')
                subplot(3,4,4);
                plot(flxO2, -zgrid, 'b', flxO2D,-zgrid,'b--',flxO2adv,-zgrid,'c--');%,flxO2D+flxO2adv,-zgrid,'r--');
                legend('tot','diff','adv','diff+adv');
                legend boxoff;
                xlabel ('O2 trsp(mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title ('O2 vert transport')
                
                
                for i=1:length(zgrid)
                    [NO3(i), flxNO3(i)] = res.zNO3.calcNO3(zgrid(i), bsd, res.swi, res);
                end
                subplot(3,4,5)
                plot(NO3, -zgrid, 'b')
                hold on
                plot([0,res.swi.NO30], [-bsd.zbio,-bsd.zbio], 'k--')
                xlabel ('NO3 (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('NO3 (mol/cm^3)')
                subplot(3,4,6)
                plot(flxNO3, -zgrid, 'b')
                xlabel ('NO3 trsp(mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title ('NO3 vert transport');
                
                
                
                for i=1:length(zgrid)
                    [NH4(i), flxNH4(i)] = res.zNH4.calcNH4(zgrid(i), bsd, res.swi, res);
                end
                subplot(3,4,7)
                plot(NH4, -zgrid, 'b')
                hold on
                plot([0,res.swi.NH40], [-bsd.zbio,-bsd.zbio], 'k--')
                hold off
                xlabel ('NH4 (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('NH4 (mol/cm^3)')
                subplot(3,4,8)
                plot(flxNH4, -zgrid, 'b');
                xlabel ('NH4 trsp(mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title ('NH4 vert transport')
                % Till here
                subplot(3,4,9)
                for i=1:length(zgrid)
                    [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), bsd, res.swi, res);
                end
                plot(SO4, -zgrid, 'b')
                hold on
                plot([0,res.swi.SO40], [-bsd.zbio,-bsd.zbio], 'k--')
                hold off
                %xlim([0 SO40])
                xlabel ('SO4 (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('SO4 (mol/cm^3)')
                subplot(3,4,10)
                plot(flxSO4, -zgrid, 'b');
                xlabel ('SO4 trsp(mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title ('SO4 vert transport')
                
                subplot(3,4,11)
                for i=1:length(zgrid)
                    [H2S(i), flxH2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                end
                plot(H2S, -zgrid, 'b')
                hold on
                plot([0,res.swi.H2S0], [-bsd.zbio,-bsd.zbio], 'k--')
                xlabel ('H2S (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('H2S (mol/cm^3)')
                subplot(3,4,12)
                plot(flxH2S, -zgrid, 'b');
                xlabel ('H2S trsp(mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title ('H2S vert transport')
                
                
                if(~swi.TwoG_OM_model)
                    figure;
                    %%% TOC
                    subplot(3,2,1)
                    for i=1:length(zgrid)
                        [C(i), C1(i,:)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                        [Cflx(i), C1flx(i,:)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                    end
                    color = parula(swi.nG);
                    for G = 1:swi.nG
                        plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
                        hold on
                    end
                    
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    
                    hold off
                    xlabel ('TOC (wt%)')
                    ylabel('Depth (cm)')
                    title('Total TOC (wt%)')
                end
            end
            
        end
        
        function plot_TOC(res, ~, swi, str_date)
            fig_toc =figure;
            
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            bsd = res.bsd;
            zgrid = 0:0.5:bsd.zinf;
            
            %% TOC
            for i=1:length(zgrid)
                [C(i), C1(i,:)] = res.zTOC_RCM.calcC( zgrid(i), bsd, res.swi, res);
                [Cflx(i), C1flx(i,:)] = res.zTOC_RCM.calcCflx( zgrid(i), bsd, res.swi, res);
            end
            subplot(121)
            color = parula(swi.nG);
            for G = 1:swi.nG
                plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
                hold on
            end
            xlabel ('TOC (wt%)')
            ylabel('Depth (cm)')
            % Plot sum (TOC)
            subplot(122)
            plot(100*C*12/bsd.rho_sed, -zgrid, ':k')
            hold on
            t=xlim;         % to draw penetration depths the correct lengths
            plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
            plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
            plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
            plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
            plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
            
            hold off
            xlabel ('TOC (wt%)')
            ylabel('Depth (cm)')
            %            title('Total TOC (wt%)')
            print(fig_toc, '-depsc2', [str_date '_TOC_PROFILES.eps']);
            
            

            
        end
        
        function plot_TOC_O2_column(res, debug, swi, str_date)
            % plot single sediment column vs depth
            
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            bsd = res.bsd;
            zgrid = 0:0.1:bsd.zinf;
            
            % CONCENTRATIONS WITHOUT PO4
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            % Nitrogen included?
            if(swi.Nitrogen)
                figure;
                
                %%% TOC
                subplot(1,2,1)
                for i=1:length(zgrid)
                    [C(i), C1(i,:)] = res.zTOC_RCM.calcC( zgrid(i), bsd, res.swi, res);
                    [Cflx(i), C1flx(i,:)] = res.zTOC_RCM.calcCflx( zgrid(i), bsd, res.swi, res);
                end
                color = parula(swi.nG);
                for G = 1:swi.nG
                    plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
                    hold on
                end
                % Plot sum (TOC)
                plot(100*C*12/bsd.rho_sed, -zgrid, ':k')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                
                hold off
                xlabel ('TOC (wt%)')
                ylabel('Depth (cm)')
                
                
                text( 3/5*(100*C(1)*12/bsd.rho_sed) ,-60,{['a:' num2str(swi.p_a)],...
                    ['SFD:' num2str(bsd.wdepth)]})
                
                %%% O2
                if(res.zox>0.0)
                    for i=1:length(zgrid)
                        [O2(i), flxO2(i), flxO2D(i), flxO2adv(i)] = res.zO2.calcO2(zgrid(i), bsd, res.swi, res);
                    end
                else
                    O2(i) = 0.0;
                end
                subplot(1,2,2)
                plot(O2, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                %                 ylim([-20 0.0])
                xlabel ('O_2 (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('O2 (mol/cm^3)')
                
            end
        end
        
        function [k, C0i, Fnonbioi] = RCM(bsd, swi)
            %% Implemented by Philip Pika
            %% Described in Pika et al. (2020) GMD
            
            % For comparison with original 2G results
            if swi.nG == 1          % use 1G to compare with Thullner ea '09
                k = swi.k;
                C0i = swi.C0_nonbio;
                Fnonbioi =  swi.C0_nonbio*(1-bsd.por)*bsd.w; % NonBioturbated SWI
                
            elseif swi.nG == 2          % nothing needs to be done here - all set in default_swi()
                %                 F(1) = 0.5;
                %                 F(2) = 0.5;
                %                 C0i = F.*swi.C0_nonbio * 1e-2/12*bsd.rho_sed;       % TOC@SWI (wt%) -> (mol/cm^3 bulk phase), 2.5 sed.density (g/cm3) 0.1
                %                 Fnonbioi = swi.C0_nonbio*(1-bsd.por)*bsd.w;        % [mol/(cm2 yr)] according non-bioturbated flux
                %                 k = [0.01 0.0001];
                %                 swi.p_a = NaN;
                %                 swi.p_nu = NaN;
            else
                emin = -15; % as in Dale ea. '15: -10;
                emax = -log10(swi.p_a)+2;  % as in Dale ea. '15: 2    % upper k limit for k-bins of multi-G approximation
                if emin >= emax;error('emin >= emax, this cannot be!');end
                %                 if emax >= log10(200);emax=log10(200);end
                
                k(1)= 10^(emin);
                kk(1)=10^(emin);
                F(1) = gammainc(swi.p_a*10^emin,swi.p_nu,'lower');
                kk(swi.nG)=10^(emax);
                k(swi.nG)=10^(emax);
                F(swi.nG) = gammainc(swi.p_a*10^emax,swi.p_nu,'upper');
                
                % Define the b.c. for all the intermediate fractions
                
                G=2:swi.nG-1;
                ne=emin+(1:swi.nG-2).*(emax-emin)./(swi.nG-1);
                kk(2:swi.nG-1)=10.^ne;
                G_inc_0 = gammainc(swi.p_a*kk(1:swi.nG-2),swi.p_nu,'upper'); % G-1 = 1:end-2
                G_inc_1 = gammainc(swi.p_a*kk(2:swi.nG-1),swi.p_nu,'upper'); % G = 2:end-1
                F(2:swi.nG-1) = (G_inc_0 - G_inc_1);
                % calculate the mean degradation rate for the intermediate OM fractions
                k(G)=kk(1:swi.nG-2)+(kk(2:swi.nG-1)-kk(1:swi.nG-2))/2;
                %                 % DH: 11.12.20 limit k-value to 1000 yr-1
                %                 k(k>1000)=1000;
                F(F<=eps)=eps;
                if abs(sum(F)-1) > 0.0001
                    warning('F~=1!!');
                end
                
                Fnonbioi = F.* ( swi.C0*(1-bsd.por)*bsd.w ); % NonBioturbated SWI
                C0i = F.*swi.C0;
%               in version before aligned with RECCAP2:                
%                if(swi.flux)
%                    Fnonbioi = F.* swi.FOM_total; % Dom was
%                else
%                  	Fnonbioi = F.*( swi.C0_nonbio*(1-bsd.por)*bsd.w ); % NonBioturbated SWI
%                end
%                C0i = F.*swi.C0_nonbio;
                
                % plot Fractions vs log(k)
                plot_fractions = false;
                if(plot_fractions)
                    benthic_test.plot_TOC_fractions(k, F)
                end
            end
        end
        
        function plot_TOC_fractions(k, F)
            
            log_k = log10(k);
            % color scale
            n_frac = size(k,2);
            color = parula(n_frac);

%             for G = 1:swi.nG
%                 plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
%             end
           
            fig_toc_fractions =figure;         

            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            box on
            hold on
            for k = 1:n_frac
            	b(k) = bar(log_k(k), F(k), 0.15, 'FaceColor',color(k,:));            
            end
            plot(log_k, F, 'k')
            xlabel ('log(k)')
            ylabel('Fraction of total OM')
            %            title('Total TOC (wt%)')
            print(fig_toc_fractions, '-depsc2', 'TOC_fractions.eps');
            
            

            
        end        
        
        function [k, C0i, Fnonbioi] = RCM_Dale_2015(bsd, swi)
            %% specify k and F speficically as in Dale ea. 2015 and use input flux to calculate non-bioturbated SWI-concentration of POC
            
            
            if(swi.Test_Dale_14G)
                
                k(1)= 1e-10;
                for i=2:13
                    j = -12+i;
                    k(i)= 3.16*10^j;
                end
                k(14)= 100;
                %                k(1:10) = 0.08;    % in case we want to degrade all
                
                F(1) = 0.021746;
                F(2) = 0.00725275;
                F(3) = 0.0096717;
                F(4) = 0.0128974;
                F(5) = 0.017199;
                F(6) = 0.0229352;
                F(7) = 0.0305846;
                F(8) = 0.0407852;
                F(9) = 0.0543879;
                F(10) = 0.0725265;
                F(11) = 0.0967046;
                F(12) = 0.12881;
                F(13) = 0.169822;
                F(14) = 0.314677;
            else
                emin = -15; % as in Dale ea. '15: -10;
                emax = -log10(swi.p_a)+2;  % as in Dale ea. '15: 2    % upper k limit for k-bins of multi-G approximation
                if emin >= emax;error('emin >= emax, this cannot be!');end
                %                 if emax >= log10(200);emax=log10(200);end
                
                k(1)= 10^(emin);
                kk(1)=10^(emin);
                F(1) = gammainc(swi.p_a*10^emin,swi.p_nu,'lower');
                kk(swi.nG)=10^(emax);
                k(swi.nG)=10^(emax);
                F(swi.nG) = gammainc(swi.p_a*10^emax,swi.p_nu,'upper');
                
                % Define the b.c. for all the intermediate fractions
                
                G=2:swi.nG-1;
                ne=emin+(1:swi.nG-2).*(emax-emin)./(swi.nG-1);
                kk(2:swi.nG-1)=10.^ne;
                G_inc_0 = gammainc(swi.p_a*kk(1:swi.nG-2),swi.p_nu,'upper'); % G-1 = 1:end-2
                G_inc_1 = gammainc(swi.p_a*kk(2:swi.nG-1),swi.p_nu,'upper'); % G = 2:end-1
                F(2:swi.nG-1) = (G_inc_0 - G_inc_1);
                k(G)=kk(1:swi.nG-2)+(kk(2:swi.nG-1)-kk(1:swi.nG-2))/2;
                F(F<=eps)=eps;
            end
            
            if abs(sum(F)-1) > 0.0001;warning('F~=1!!');end
            
            % use the input flux of 6.2 mmol m-2 d-1 from Dale
            % BC Dale ea 2015:
            %                POC_flux = [0.5 1 2 4 6 8 10 12 14 16];
            
            FPOC_Dale = swi.POC_flux(swi.POCi)*10^(-3)/100^2*365;   %   in mol / (cm^2*yr)
            %                FPOC_Dale = 6.0*10^(-3)/100^2*365;   %   in mol / (cm^2*yr)
            Fnonbioi = F.* FPOC_Dale; % NonBioturbated SWI
            C0i = Fnonbioi/((1-bsd.por)*bsd.w); % this value is not used anywhere - C0i is calculated from SWI-flux in benthic_zTOC_RCM.m
            %                Fnonbioi = F.* ( swi.C0_nonbio*(1-bsd.por)*bsd.w ); % NonBioturbated SWI
            %                C0i = F.*swi.C0_nonbio;
            
        end
        
        function run_PO4flux_SA_OMEN()
            % make PO4-SWI flux SA for changing boundary conditions in Corg and O2
            clear
            
            swi=benthic_test.default_swi()
            %            % set date-time
            %            str_date = datestr(now,'ddmmyy_HH_MM_SS');
            for i=1:11   % 51
                swi.C01 = (0.1+(i-1)*0.2)*1e-2/12*2.5;     % 2.5 is rho_sed                      % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
                swi.C02 = swi.C01;
                Corg(i)=2*(0.1+(i-1)*0.2);
                for j=1:21  % 51
                    swi.O20=(j-1)*3.0E-009;                                 % O2  concentration at SWI (mol/cm^3)
                    res=benthic_test.test_benthic(1,swi);
                    SWI_PO4(i,j)=res.flxswi_P;
                    O20(j)=(j-1)*3.0E-009;
                    %                  if((i==2 && j==2) || (i==3 && j==4))2.78E-005
                    
                    %                      benthic_test.plot_column(res, false, swi, '00_expl_profiles')
                    %                  end
                end
            end
            
            figure;
            hold on
            [C,h] = contourf(O20,Corg,SWI_PO4);
            clabel(C,h,'FontSize',16);
            colorbar();
            box on
            hold off
            xlabel('Ocean O_2')
            ylabel ({'Corg'}); %;'(\mumol cm^{-2}yr^{-1})'})
            print('-depsc', '0_PO4_SWI-flux_SA_1811_NoAnoxicChange_FullP-cycle_higherPO4a_FasterDegrad_1000m');
            %            print('-depsc', 'PO4_SWI-flux_SA_0907_4PO4_006500001_15Corg_150O2_3000m_ksfast');
            
        end
        
        
    end
    
end
