%% Make sure every grid-cell with TOC-values also has a value for the BC
%% ALSO create the total FeOOH flux from the aligned sedimentation rate
% here using the a resolution of 1 degrees
% use mean of neighboring cells -
% might need to run it multiple times to
% assure all grid-cells are filled with a value
% at the end make grid-cells NaN that are NaN in Lee-TOC
clear all

path(path,'/home/domhu/Documents/MATLAB/M_Map');

align_Bcs = true;
plot_original_BC = false;    % also plot the original BC without filled values
make_feooh_from_sedrate = true;

str_date = [datestr(date,7), datestr(date,5), datestr(date,11)];

% load TOC data and lat and long
load('./BC_1degree/Lee_toc_lr_weighted.mat')
toc = Lee_toc_lr_weighted;
load('./BC_1degree/long_lr.mat')
load('./BC_1degree/lat_lr.mat')

% load BC in 1degree resolution
load './data/WOA2018/Tmp_BW_2018_01.mat'
Tmp_BW_WOA2018_lr = Tmp_BW_2018_01;
load './data/WOA2018/O2_BW_WOA2018.mat'
O2_BW_WOA2018_lr = O2_BW_WOA2018;
load './data/WOA2018/NO3_BW_2018.mat'
NO3_BW_WOA2018_lr = NO3_BW_2018;

load('./BC_1degree/sed_lr_weighted.mat');

% seafloor depth from NASA: https://neo.sci.gsfc.nasa.gov/view.php?datasetId=GEBCO_BATHY
load './data/GEBCO_BATHY_2002_360x180.CSV'
water_depth = flipud(GEBCO_BATHY_2002_360x180);
water_depth_nan=find(water_depth>=0);
water_depth(water_depth_nan)=NaN;

% todo
% load 'feooh_Thullner_TOTAL_BW.mat'
% load 'FluxFe_total_Burw_BW.mat'
% load 'FluxFe_total_Burw_Daleunits_BW.mat'

toc_notNaN = (~isnan(toc));
size_toc = sum(sum(toc_notNaN));
[n,m]=size(toc);


% loop through the different boundary consitions (i.e. 5 different ones)
if(align_Bcs)
    for bc=5:5
        switch bc
            case 1  % Temperature
                BCond_BW_lr_ori = Tmp_BW_WOA2018_lr;
                BCond_BW_lr = Tmp_BW_WOA2018_lr;
            case 2  % O2
                BCond_BW_lr_ori = O2_BW_WOA2018_lr;
                BCond_BW_lr = O2_BW_WOA2018_lr;
            case 3  % NO3
                BCond_BW_lr_ori = NO3_BW_WOA2018_lr;
                BCond_BW_lr = NO3_BW_WOA2018_lr;
            case 4  % water-depth
                BCond_BW_lr_ori = water_depth;
                BCond_BW_lr = water_depth;
            case 5  % Holocene sedimentation rate (cm/yr)
                BCond_BW_lr_ori = sed_lr_weighted;
                BCond_BW_lr = sed_lr_weighted;
        end
        BCond_notNaN = (~isnan(BCond_BW_lr));
        size_BCond_orig = sum(sum(BCond_notNaN));
        size_BCond_new = size_BCond_orig;   % inituialize for while-loop
        size_BCond_old = 0;
        iter = 0;
        
        while size_BCond_new ~= size_BCond_old
            iter=iter+1
            % --- could improve the loop-through-all by only going to the
            % cells where isnan(BC) but a value for TOC exists
            for i=2:n-1         % what about first & last line/row? -> handle separate!?
                for j=1:m
                    if(toc_notNaN(i,j))
                        if(~BCond_notNaN(i,j))    % if NaN calc mean of neighbors
                            if(j==1)        % first column - use last column to calculate mean
                                BCond_BW_lr(i,j)=mean([BCond_BW_lr(i+1,j),BCond_BW_lr(i,j+1),BCond_BW_lr(i-1,j),BCond_BW_lr(i,m)],'omitnan');
                            elseif(j==m)    % last column - use first column to calculate mean
                                BCond_BW_lr(i,j)=mean([BCond_BW_lr(i+1,j),BCond_BW_lr(i,1),BCond_BW_lr(i-1,j),BCond_BW_lr(i,j-1)],'omitnan');
                            else
                                BCond_BW_lr(i,j)=mean([BCond_BW_lr(i+1,j),BCond_BW_lr(i,j+1),BCond_BW_lr(i-1,j),BCond_BW_lr(i,j-1)],'omitnan');
                            end
                        end
                    end
                    
                end
            end
            size_BCond_old=size_BCond_new;
            BCond_notNaN_new = (~isnan(BCond_BW_lr));
            size_BCond_new = sum(sum(BCond_notNaN_new));
        end
        
        % check that all TOC grid-cells have a value in the BC files:
        % and plot later where values in BC are missing
        where_toc_notNaN = find(~isnan(toc));
        BCond_help = BCond_BW_lr(where_toc_notNaN);
        BCond_nan = (isnan(BCond_help));
        size_BCond_NaN(bc,:) = size(BCond_help(BCond_nan));     % save how many cells with a TOC-value are missing this BC
        
        % check where the data is missing:
        BCond_isNaN = (isnan(BCond_BW_lr));     % isnan in BC
        Loc_BC_missing = toc_notNaN .* BCond_isNaN;
        
        [Loc_BC_miss_row, Loc_BC_miss_colum] = find(Loc_BC_missing==1);
        % get the coordinates for missing missing BC
        Loc_BC_miss_lat = lat_lr(Loc_BC_miss_row);
        Loc_BC_miss_long = long_lr(Loc_BC_miss_colum);
        
        switch bc
            case 1
                Tmp_BW_lr_aligned = BCond_BW_lr;
                save('./BC_1degree/Tmp_BW_lr_aligned.mat' , 'Tmp_BW_lr_aligned')
                txt_file = 'Temp_WOA_';
                txt_title = 'Temperature at the SWI (\circC)';
                levels = [-3:1:31];
            case 2  % O2
                O2_BW_WOA2018_lr_aligned = BCond_BW_lr;
                save('./BC_1degree/O2_BW_WOA2018_lr_aligned.mat' , 'O2_BW_WOA2018_lr_aligned')
                txt_file = 'O2_muM_WOA2018_';
                txt_title = 'O_2 at the SWI (\muM)';
                levels = [0:10:400];
            case 3  % NO3
                NO3_BW_lr_aligned = BCond_BW_lr;
                save('./BC_1degree/NO3_BW_lr_aligned.mat' , 'NO3_BW_lr_aligned')
                txt_file = 'NO3_muM_WOA_';
                txt_title = 'NO_3 at the SWI (\muM)';
                levels = [0:1:52];
            case 4  % water-depth
                water_depth_aligned = BCond_BW_lr;
                save('./BC_1degree/water_depth_aligned.mat' , 'water_depth_aligned')
                txt_file = 'Depth_NASA_';
                txt_title = 'Seafloor depth (m)';
                levels=[-8000:200:0];
            case 5  % Sedimentation rate (cm/yr) -- one could also use the Burwicz parameterization for cells that are missing instead of the mean
                sed_lr_weighted_aligned = BCond_BW_lr;
                save('./BC_1degree/sed_lr_weighted_aligned.mat' , 'sed_lr_weighted_aligned')
                txt_file = 'sed_lr_weighted';
                txt_title = 'Sedimentation rate (cm/yr)';
                levels = [0:0.005:0.13];
        end
        
        % Robinson projection
        fig_SWI_BCond(bc) = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        % set(gca,'FontSize',30)
        %        depth_levels=[-8000:200:0];
        [C,h] = m_contourf(long_lr, lat_lr, BCond_BW_lr, levels);  %,depth_levels);
        set(h,'LineColor','none')
        title(txt_title)
        m_coast('linewidth',0.5); %('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        % plot missing temp-grids
        %        scatter(45, 45, 100, 'xr')
        %        scatter(Loc_BC_miss_long, Loc_BC_miss_lat, 100, 'xr')
        m_plot(Loc_BC_miss_long, Loc_BC_miss_lat,'xr','MarkerSize',8);
        m_text(45, 60, ['missing ' int2str(size_BCond_NaN(bc,1)), ' of ', int2str(size_toc)],'FontSize',6);
        % shading interp;
        colorbar ('horizontal')
        colormap(parula)
        %       caxis([-5000.0 0.0])
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_SWI_BCond(bc),'-depsc2', ['./SWI_eps/SWI_' txt_file, str_date '_lr_1deg_FILLED.eps']);
        
        if(plot_original_BC)
            % Robinson projection
            fig_SWI_BCond_ori(bc) = figure;
            m_proj('Robinson','longitudes',[-180 179.99], ...
                'latitudes',[-90 90]);
            hold on;
            [C,h] = m_contourf(long_lr, lat_lr, BCond_BW_lr_ori, levels);  %,depth_levels);
            set(h,'LineColor','none')
            title([ txt_title , ' original'])
            m_coast('linewidth',0.5); %('patch',[.7 .7 .7]);
            m_grid('fontsize',8);
            % hold on;
            % shading interp;
            colorbar ('horizontal')
            colormap(parula)
            %       caxis([-5000.0 0.0])
            xlabel('Longitude')
            ylabel('Latitude')
            print(fig_SWI_BCond_ori(bc),'-depsc2', ['./SWI_eps/SWI_' txt_file, str_date '_lr_1deg_FILLED_ori.eps']);
        end
        
        
    end
end


if(make_feooh_from_sedrate)
    % Use the aligned sedimentation rate from above to
    % Calculate the total Fe flux as fct. of Holocene sedimentation rate
    % compare Dale et al. (2015) Table 2:
    % FluxFe_total = frac_fe * w * (1-por)*rho/Aw   % in mol /(cm2 yr)
    
    %    load sed_holo.mat   % in cm/yr  -- 1/4 degree resolution
    % load the updated sedimentation rate:
    load('./BC_1degree/sed_lr_weighted_aligned.mat');
    
    w = sed_lr_weighted_aligned;  % in cm/yr
    frac_fe = 0.05;     %  Fe content in average sedimentary rock (~5%) [Garrels and Mackenzie, 1971] which is similar to Fe content in red clays [Glasby, 2006].
    por = 0.85;          % porosity of compacted sediment
    rho = 2.5;          % dry sediment density, in g cm-3)
    Aw = 55.8;          %  standard atomic weight of iron (55.8 g mol-1)
    
    FluxFe_total_Burw_BW_01_por85_aligned = frac_fe * w * (1-por)*rho/Aw *10^6 ;                %  in mumol/(cm2 yr)
    FluxFe_total_Burw_Daleunits_BW_01_por85_aligned = FluxFe_total_Burw_BW_01_por85_aligned *100^2 /365;    % in mumol/(m2 d))
    
    por = 0.74;          % avg value as in Wallmann et al. (2012)
    FluxFe_total_Burw_BW_01_por75 = frac_fe * w * (1-por)*rho/Aw *10^6 ;                %  in mumol/(cm2 yr)
    FluxFe_total_Burw_Daleunits_BW_01_por74 = FluxFe_total_Burw_BW_01_por75 *100^2 /365;    % in mumol/(m2 d))
    
    % save the BW concentrations
    save('./BC_1degree/FluxFe_total_Burw_BW_aligned.mat' , 'FluxFe_total_Burw_BW_01_por85_aligned' , 'FluxFe_total_Burw_BW_01_por75')
    save('./BC_1degree/FluxFe_total_Burw_Daleunits_BW_aligned.mat' , 'FluxFe_total_Burw_Daleunits_BW_01_por85_aligned', 'FluxFe_total_Burw_Daleunits_BW_01_por74')
    
    %     % 1/4 degree resolution
    %     load 'lat.dat'
    %     load 'long.dat'
    
    % 1 degree resolution
    load('./BC_1degree/lat_lr.mat');
    load('./BC_1degree/long_lr.mat');
    lat = lat_lr;
    long = long_lr;
    
    % Robinson projection
    fig_feooh = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    levels=[0:100:2200];
    [C,h] = m_contourf(long, lat, FluxFe_total_Burw_Daleunits_BW_01_por85_aligned,levels);
    set(h,'LineColor','none')
    title('Total flux of Fe(OH)_3 (\mumol m-2 d-1)')
    m_coast('linewidth',0.5); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 1300.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_feooh,'-depsc2', ['FeOOH_flux_por085_Daleunits.eps']);
    
    % Robinson projection
    fig_feooh_por74 = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    levels=[0:100:2200];
    [C,h] = m_contourf(long, lat, FluxFe_total_Burw_Daleunits_BW_01_por74,levels);
    set(h,'LineColor','none')
    title('Total flux of Fe(OH)_3 (\mumol m-2 d-1)')
    m_coast('linewidth',0.5); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 1300.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_feooh_por74,'-depsc2', ['FeOOH_flux_por074_Daleunits.eps']);
    
    
    
    
end
