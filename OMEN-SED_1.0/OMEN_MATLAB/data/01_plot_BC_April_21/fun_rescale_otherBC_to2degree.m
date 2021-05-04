%% Reduce the 1 degree resolution files (that have been aligned with Lee-TOC) created earlier to 2 degrees
% assure all grid-cells with TOC are filled with a value -- probably don't need to
% align as 1-degree data has been aligned with Lee-TOC

clear all

path(path,'/home/domhu/Documents/MATLAB/M_Map');

reduce_all_BCs = true;  % reduce other BC from 1 to 2 degrees
align_Bcs = false;              % probably not needed because 1 degree has been aligned
make_param_a = false;   % at the send from updated sed-holo

str_date = [datestr(date,7), datestr(date,5), datestr(date,11)];

% O2_filename='./WOA2018/woa18_all_o00_01.nc';
% long_WOA_01 = ncread(O2_filename,'lon');
% lat_WOA_01 = ncread(O2_filename,'lat');

if(reduce_all_BCs)
    % load Lee-TOC in 2-degrees
    load './BC_2degree/Lee_toc_lr_weighted.mat';
    
    %% load data in 1 degree resolution    
    load './BC_1degree/Tmp_BW_lr_aligned.mat'                             % BW temperature [°C]
    BC{1} = Tmp_BW_lr_aligned;
    load './BC_1degree/sed_lr_weighted_aligned.mat'                           % Sedimentation rate [cm/yr] -- as in Bradley ea. 2020, after Burwicz ea. 2011
    BC{2} = sed_lr_weighted_aligned;
    load './BC_1degree/O2_BW_WOA2018_lr_aligned.mat'                  	% BW  O2 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
    BC{3} = O2_BW_WOA2018_lr_aligned;
    load './BC_1degree/NO3_BW_lr_aligned.mat'                             % BW  NO3 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
    BC{4} = NO3_BW_lr_aligned;
    load './BC_1degree/FluxFe_total_Burw_BW_aligned.mat'               % Total FeOOH settling-flux [mumol/(cm2 yr)] -- just 16.67% of this is available for OM-degradation; 
    BC{5} = FluxFe_total_Burw_BW_01_por85_aligned;                      % used por=0.85 to calculate it
    BC{6} = FluxFe_total_Burw_BW_01_por75;                              % used por=0.75 to calculate it
    load './BC_1degree/water_depth_aligned.mat'                           % Seafloor depth [m]
    BC{7} = water_depth_aligned;
    load './BC_1degree/a_value_total_PK.mat';                           % Philip a-values [yr-1]
    BC{8} = Parameter_a_2Dtotal;
    load './BC_1degree/a_lr_aligned_1to100.mat';                           % a-value using Sandra's parameterization rescales to 1 - 100[yr-1]
    BC{9} = a_lr_aligned;
    %                 load ./BC_1degrees/long_WOA_01.mat
    %                 long_01_in = long_WOA_01;
    %                 load ./BC_1degrees/lat_WOA_01.mat
    %                 lat_01_in = lat_WOA_01;
    
    %% load 2 degree lat/long
    load './BC_2degree/lat_lr.mat'
    load './BC_2degree/long_lr.mat'
    
    
    [n,m]=size(Tmp_BW_lr_aligned);
    
    for k=9:9 %size(BC,2)
        BC_in = BC{k};
        
        x_new = 0;
        y_new = 0;
        
        for x=1:2:n-1
            x_new=x_new+1;
            y_new = 0;
            for  y=1:2:m-1
                y_new = y_new+1;
                % BC
                BC_sub = BC_in(x:x+1, y:y+1);               % 2x2 submatrix that will be combined
                BC_nanmean =  nanmean(nanmean(BC_sub));     % calculate mean of toc
                BC_02(x_new, y_new) = BC_nanmean;            % save in new matrix
            end
        end
        
        % write NaN where Lee-TOC is NaN:
        TOC_isNaN = isnan(Lee_toc_lr_weighted);
        BC_02(TOC_isNaN)= NaN;
        
        BC_02_out{k} = BC_02;
        
        switch k
            case 1
                Tmp_BW_02_updated = BC_02;
                save('./BC_2degree/Tmp_BW_02_updated.mat' , 'Tmp_BW_02_updated')
                txt_file = 'Temp_WOA_02_';
                txt_title = 'Temperature at the SWI (\circC) -- 2^\circ';
                levels = [-5:0.5:25];
                limits = [0 25];
            case 2
                sed_holo_02 = BC_02;
                save('./BC_2degree/sed_holo_02.mat' , 'sed_holo_02')
                txt_file = 'Sed_holo_02_';
                txt_title = 'Sedimentation rate (cm/yr) -- 2^\circ';
                levels = [0:0.005:0.14];
                limits = [0 0.12];
            case 3  % O2
                O2_BW_WOA2018_02_updated = BC_02;
                save('./BC_2degree/O2_BW_WOA2018_02_updated.mat' , 'O2_BW_WOA2018_02_updated')
                txt_file = 'O2_muM_WOA2018_02_';
                txt_title = 'O_2 at the SWI (\muM) -- 2^\circ';
                levels = [0:1.0:500];
                limits = [0 350];
            case 4  % NO3
                NO3_BW_02_updated = BC_02;
                save('./BC_2degree/NO3_BW_02_updated.mat' , 'NO3_BW_02_updated')
                txt_file = 'NO3_muM_WOA_02_';
                txt_title = 'NO_3 at the SWI (\muM) -- 2^\circ';
                levels = [0:0.5:80];
                limits = [0 50];
            case 5  % Total feooh 0.85 - after Burwicz in mumol/(cm2 yr)
                FluxFe_total_por085_Burw_BW_02_updated = BC_02;                
                save('./BC_2degree/FluxFe_total_por085_Burw_BW_02_updated.mat' , 'FluxFe_total_por085_Burw_BW_02_updated')
                txt_file = 'FluxFe_total_por085_Burw_BW_02_updated';
                txt_title = 'Total flux of Fe(OH)_3 por = 0.85 (\mumol/(cm^2 yr)) -- 2^\circ';
                levels = [0:0.5:45];
                limits = [0 30];
            case 6  % Total feooh 0.75 - after Burwicz in mumol/(m2 d))
                FluxFe_total_por075_Burw_BW_02_updated = BC_02;
                save('./BC_2degree/FluxFe_total_por075_Burw_BW_02_updated.mat' , 'FluxFe_total_por075_Burw_BW_02_updated')
                txt_file = 'FluxFe_total_por075_Burw_BW_02_updated_';
                txt_title = 'Total flux of Fe(OH)_3 por=0.75 (\mumol/(cm^2 yr)) -- 2^\circ';
                levels = [0:0.5:45];
                limits = [0 30];
            case 7  % water-depth
                water_depth_02_updated = BC_02;
                save('./BC_2degree/water_depth_02_updated.mat' , 'water_depth_02_updated')
                txt_file = 'Depth_NASA_02_';
                txt_title = 'Seafloor depth (m) -- 2^\circ';
                levels = [-8000:100:0];
                limits = [-5000 0];
            case 8  % Philip's a-values
                a_values_PK_updated = BC_02;
                save('./BC_2degree/a_values_PK_updated.mat' , 'a_values_PK_updated')
                txt_file = 'A_values_PK_';
                txt_title = 'a-value (yr) -- 2^\circ';
                levels = [0:1:100];
                limits = [0 50];
            case 9  % a-value using Sandra's parameterization
                a_values_SA_1to100_updated = BC_02;
                save('./BC_2degree/a_values_SA_1to100_updated.mat' , 'a_values_SA_1to100_updated')
                txt_file = 'A_values_1to100_';
                txt_title = 'a-value (yr) -- 2^\circ';
                levels = [0:1:100];
                limits = [0 100];
        end
        
        %% plot BC
        % using Robinson projection
        fig_SWI_BCond(k) = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
%        sed_levels = [0:0.005:0.13];
        [C,h] = m_contourf(long_lr, lat_lr, BC_02, levels);  
        set(h,'LineColor','none')
        title(txt_title)
        m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_SWI_BCond(k),'-depsc2', ['./BC_2degree/SWI_' txt_file, str_date '_02.eps']);
        
    end
                        
    % plot Lee-toc in 2 degrees
        % using Robinson projection
        fig_LeeTOC = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        sed_levels = [0:0.1:5];
        [C,h] = m_contourf(long_lr, lat_lr, Lee_toc_lr_weighted, sed_levels); 
        set(h,'LineColor','none')
        title('Lee TOC at SWI (wt%) -- 2^\circ')
        m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis([0.0 3.0])
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_LeeTOC,'-depsc2', ['./BC_2degree/SWI_LeeTOC_' str_date '_02.eps']);
    
    
end
