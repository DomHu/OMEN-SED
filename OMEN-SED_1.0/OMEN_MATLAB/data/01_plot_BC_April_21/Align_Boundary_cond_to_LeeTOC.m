%% Make sure every grid-cell with TOC-values also has a value for the BC
% here using the a resolution of 1/4 degrees
% use mean of neighboring cells -
% (Note: For sedimentation rate use the Burwicz et al. (2011) equation)
% might need to run it multiple times to
% assure all grid-cells are filled with a value
% at the end make grid-cells NaN that are NaN in Lee-TOC
clear all

align_Bcs = true;
make_param_a = false;
str_date = [datestr(date,7), datestr(date,5), datestr(date,11)];

% load TOC data and lat and long
load('./BC_Quarterdegree/Lee_toc_lr_weighted.mat')
toc = Lee_toc_lr_weighted;
load('./BC_Quarterdegree/long_lr.mat')
load('./BC_Quarterdegree/lat_lr.mat')

% load BC in 1/4 resolution
load './data/WOA2018/Tmp_BW_WOA2018_hr.mat'
load './data//WOA2018/O2_BW_WOA2018_hr.mat'
load './data/WOA2018/NO3_BW_WOA2018_hr.mat'

load('./BC_Quarterdegree/sed_holo_clean.mat');  % this was provided in 1/4 resolution

% todo
% load 'water_depth.mat'
% load 'feooh_Thullner_TOTAL_BW.mat'
% load 'FluxFe_total_Burw_BW.mat'
% load 'FluxFe_total_Burw_Daleunits_BW.mat'

toc_notNaN = (~isnan(toc));
size_toc = sum(sum(toc_notNaN));
[n,m]=size(toc);


% loop through the different boundary consitions (i.e. 8 different ones)
if(align_Bcs)
    for bc=1:1
        switch bc
            case 1  % Temperature
                BCond_BW_hr = Tmp_BW_WOA2018_hr;
            case 2  % O2
                BCond_BW_hr = O2_BW_WOA2018_hr;
            case 3  % NO3
                BCond_BW_hr = NO3_BW_WOA2018_hr;
            case 4  % water-depth
                BCond_BW_hr = water_depth;
            case 5  % Total feooh - from Thullner
                BCond_BW_hr = feooh_Thullner_TOTAL_BW;
            case 6  % Total feooh - after Burwicz in mumol/(cm2 yr)
                BCond_BW_hr = FluxFe_total_Burw_BW;
            case 7  % Total feooh - after Burwicz in mumol/(m2 d))
                BCond_BW_hr = FluxFe_total_Burw_Daleunits_BW;
            case 8  % Holocene sedimentation rate (cm/yr)
                BCond_BW_hr = sed_holo;
        end
        BCond_notNaN = (~isnan(BCond_BW_hr));
        size_BCond_orig = sum(sum(BCond_notNaN));
        size_BCond_new = size_BCond_orig;   % inituialize for while-loop
        size_BCond_old = 0;
        iter = 0;
        
        while size_BCond_new ~= size_BCond_old
            iter=iter+1
            % --- could improve the loop-through-all by only going to the
            % cells where isnan(BC) but a value for TOC exists
            for i=2:n-1         % what about first & last line/row? -> handle separate!?   
                for j=2:m-1
                    if(toc_notNaN(i,j))
                        if(~BCond_notNaN(i,j))    % if NaN calc mean of neighbors
                            BCond_BW_hr(i,j)=mean([BCond_BW_hr(i+1,j),BCond_BW_hr(i,j+1),BCond_BW_hr(i-1,j),BCond_BW_hr(i,j-1)],'omitnan');
                        end
                    end
                    
                end
            end
            size_BCond_old=size_BCond_new;
            BCond_notNaN_new = (~isnan(BCond_BW_hr));
            size_BCond_new = sum(sum(BCond_notNaN_new));
        end
        
        % check that all TOC grid-cells have a value in the BC files:
        
        where_toc_notNaN = find(~isnan(toc));
        BCond_help = BCond_BW_hr(where_toc_notNaN);
        BCond_nan = (isnan(BCond_help));
        size_BCond_NaN(bc,:) = size(BCond_help(BCond_nan));
        
        
        switch bc
            case 1
                Tmp_BW_hr_updated = BCond_BW_hr;
                save('Tmp_BW_hr_updated.mat' , 'Tmp_BW_hr_updated')
                txt_file = 'Temp_WOA_';
                txt_title = 'Temperature at the SWI (\circC)';
            case 2  % O2
                O2_BW_WOA2018_hr_updated = BCond_BW_hr;
                save('O2_BW_WOA2018_hr_updated.mat' , 'O2_BW_WOA2018_hr_updated')
                txt_file = 'O2_muM_WOA2018_';
                txt_title = 'O_2 at the SWI (\muM)';
            case 3  % NO3
                NO3_BW_hr_updated = BCond_BW_hr;
                save('NO3_BW_hr_updated.mat' , 'NO3_BW_hr_updated')
                txt_file = 'NO3_muM_WOA_';
                txt_title = 'NO_3 at the SWI (\muM)';
            case 4  % water-depth
                water_depth_updated = BCond_BW_hr;
                save('water_depth_updated.mat' , 'water_depth_updated')
                txt_file = 'Depth_NASA_';
                txt_title = 'Seafloor depth (m)';
            case 5  % feooh Thullner
                feooh_Thullner_TOTAL_BW_hr_updated = BCond_BW_hr;
                save('feooh_Thullner_TOTAL_BW_hr_updated.mat' , 'feooh_Thullner_TOTAL_BW_hr_updated')
                txt_file = 'FeOOH_Thullner_Total_';
                txt_title = 'Thullner Total flux of Fe(OH)_3 (\mumol/(cm^2 yr))';
            case 6  % Total feooh - after Burwicz in mumol/(cm2 yr)
                FluxFe_total_Burw_BW_hr_updated = BCond_BW_hr;
                save('FluxFe_total_Burw_BW_hr_updated.mat' , 'FluxFe_total_Burw_BW_hr_updated')
                txt_file = 'FeOOH_Burwicz_Total_';
                txt_title = 'Total flux of Fe(OH)_3 (\mumol/(cm^2 yr))';
            case 7  % Total feooh - after Burwicz in mumol/(m2 d))
                FluxFe_total_Burw_Daleunits_BW_hr_updated = BCond_BW_hr;
                save('FluxFe_total_Burw_Daleunits_BW_hr_updated.mat' , 'FluxFe_total_Burw_Daleunits_BW_hr_updated')
                txt_file = 'FeOOH_Burwicz_Total_Daleunits';
                txt_title = 'Total flux of Fe(OH)_3 (\mumol/(m^2 d))';
            case 8  % Holocene sedimentation rate (cm/yr)
                sed_holo_hr_updated = BCond_BW_hr;
                save('sed_holo_hr_updated.mat' , 'sed_holo_hr_updated')
                txt_file = 'sed_holo_hr_updated';
                txt_title = 'Holocene sedimentation rate (cm/yr)';
        end
        
        fig_SWI_BCond(bc) = figure;
        set(gca,'FontSize',30)
        pcolor(long, lat, BCond_BW_hr);
        title(txt_title)
        hold on;
        shading interp;
        %contour(long, lat, toc,'LineColor','k')
        colorbar ('horizontal')
        colormap(parula)
        %caxis([-5000 0.0])
        
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_SWI_BCond(bc),'-depsc2', ['./SWI_eps/SWI_' txt_file, str_date '_HR_FILLED.eps']);
        
    end
end

if(make_param_a)
    % calculate a as fct of sedimentation rate dependent formulation for a (Arndt et al., 2013)
    
    load sed_holo_hr_updated.mat       % in cm/yr
    
    loga=3.35-14.81*sed_holo_hr_updated;                          %sedimentation rate dependent formulation for a (Arndt et al., 2013)
    a_hr_updated=10.^(loga);
    
    save('a_hr_updated.mat' , 'a_hr_updated')
    
    txt = 'log(a)=3.35-14.81\cdot w';
    
    fig_SWI_a = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, a_hr_updated);
    title('Parameter a at the SWI (yr)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    %    caxis([0.0 5.0])
    text(20,-80,txt);
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_SWI_a,'-depsc2', ['./SWI_eps/SWI_parameter_a_' str_date '_HR_FILLED.eps']);
    
    
end