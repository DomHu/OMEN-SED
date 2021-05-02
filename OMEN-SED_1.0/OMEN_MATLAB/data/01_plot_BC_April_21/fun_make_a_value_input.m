%% Make a global map for a-values
% First just use Sandra's sed-accumulation rate related parameterization
% scaled to different ranges (set them in line 43, 44)
% Second use Philip's global a-value map:
% Idea for what to do with a-values from Philip:
% Check which grid-cells are missing?
% For deeper ocean use mean of neighboring cells
% If values are consistently missing for depths shallower than x meters
% use a lower a-value -- check literature and check that final total OM
% degradation rate is reasonable
path(path,'/home/domhu/Documents/MATLAB/M_Map');

str_date = [datestr(date,7), datestr(date,5), datestr(date,11)];

clear all;

make_param_a_from_sedrate = true;
plot_fig = false;

% load data in 1 degree resolution
load('./BC_1degree/lat_lr.mat');
load('./BC_1degree/long_lr.mat');
lat = lat_lr;
long = long_lr;
load('./BC_1degree/water_depth_aligned.mat');

load('./BC_1degree/Lee_toc_lr_weighted.mat')
toc = Lee_toc_lr_weighted;



if(make_param_a_from_sedrate)
    
    % calculate a as fct of sedimentation rate dependent formulation for a (Arndt et al., 2013)
    
    load('./BC_1degree/sed_lr_weighted_aligned.mat');       % in cm/yr
    
    loga=3.35-14.81*sed_lr_weighted_aligned;                          %sedimentation rate dependent formulation for a (Arndt et al., 2013)
    a_lr_aligned_ori=10.^(loga);
    
    a_max = max(max(a_lr_aligned_ori));
    a_min = min(min(a_lr_aligned_ori));
    
    % Set new range of a-values:
    a_max_new = 100;
    a_min_new = 10.0;   % 0.5 does not seem to work; high sed-rate cells have no value
    a_lr_aligned = (a_lr_aligned_ori - a_min)./(a_max-a_min).*(a_max_new - a_min_new) + a_min_new;
    
    loga_rescaled = log10(a_lr_aligned);
    
    save(['BC_1degree/a_lr_aligned_', num2str(a_min_new) , 'to', num2str(a_max_new) ,'.mat'] , 'a_lr_aligned')
    
    txt = 'log(a)=3.35-14.81\cdot w';
    

    fig_SWI_a = figure;
m_proj('Robinson','longitudes',[-180 179.99], ...
           'latitudes',[-90 90]);
hold on; 
% set(gca,'FontSize',30)
% plot log(a)
%a_levels = [-0.1:0.01:2.0];
%[C,h] = m_contourf(long, lat, loga_rescaled, a_levels);
%plot a:
a_levels = [1:0.25:100.0];
[C,h] = m_contourf(long, lat, a_lr_aligned, a_levels);
set(h,'LineColor','none')
% title(['log(a) - ', num2str(a_min_new) , ' to ', int2str(a_max_new)])
title(['a-value: ', num2str(a_min_new) , ' to ', int2str(a_max_new)])
m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]); 
m_grid('fontsize',8);
m_text(45, 60, txt,'FontSize',6);
% hold on;
% shading interp; 
colorbar ('horizontal')
colormap(parula)
caxis([0.0 100.0])
xlabel('Longitude')
ylabel('Latitude')
print(fig_SWI_a,'-depsc2', ['SWI_eps/SWI_loga_', num2str(a_min_new) , 'to', num2str(a_max_new) ,'.eps']);


end


if(plot_fig)
    load('./data/Pika/Parameter_log10_a_2D.mat');
    Parameter_a_2D = flipud(Parameter_a_2D);
    Parameter_a_2Dtotal = 10.^Parameter_a_2D;
    
    fig01 = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    a_levels = [-0.1:0.1:2.0];
    [C,h] = m_contourf(long, lat, Parameter_a_2D, a_levels);
    set(h,'LineColor','none')
    %m_contourf(long_WOA_01, lat_WOA_01, Parameter_a_2Dtotal);
    title('log(a)')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    %caxis([0.0 25.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig01,'-depsc2', ['SWI_a_values_PK_orig.eps']);
    
    fig02 = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    levels=[-8000:200:0];
    [C,h] = m_contourf(long, lat, water_depth_aligned, levels);
    set(h,'LineColor','none')
    %m_contourf(long_WOA_01, lat_WOA_01, Parameter_a_2Dtotal);
    title('water depth')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    %caxis([0.0 25.0])
    xlabel('Longitude')
    ylabel('Latitude')
    
    % Robinson projection
    fig_toc_Lee = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    toc_levels = [0:0.1:5];
    [C,h] = m_contourf(long, lat, toc,toc_levels);
    set(h,'LineColor','none')
    title('Lee SWI-TOC (wt%)')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 3.0])
    xlabel('Longitude')
    ylabel('Latitude')
    
    
    load('./BC_1degree/sed_lr_weighted_aligned.mat');
    % using Robinson projection
    fig_sedrate_Brad = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    sed_levels = [0:0.005:0.13];
    [C,h] = m_contourf(long, lat, sed_lr_weighted_aligned, sed_levels);    % 1/4 degree resolution
    set(h,'LineColor','none')
    title('Sedimentation rate as Bradley (cm/yr)')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 0.12])
    xlabel('Longitude')
    ylabel('Latitude')
    
    
    load('data/Pika/Province_Coordinate_0.25deg.mat')
    Provinces = flipud(Provinces);
    Provinces(Provinces==0) = NaN;
    load('BC_Quarterdegree/lat_lr.mat');
    load('BC_Quarterdegree/long_lr.mat');
    
    [Loc_SouthPac_row, Loc_SouthPAc_colum] = find(Provinces==4);
    % get the coordinates for missing missing BC
    Loc_SouthPac_lat = lat_lr(Loc_SouthPac_row);
    Loc_SouthPac_long = long_lr(Loc_SouthPAc_colum);
    
    
    
    % Robinson projection
    fig_toc_provinces = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    provinces_levels = [0:1:22];
    [C,h] = m_contourf(long_lr, lat_lr,Provinces,provinces_levels);
    set(h,'LineColor','none')
    title('Provinces')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    m_plot(Loc_SouthPac_long, Loc_SouthPac_lat,'xr','MarkerSize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    m_colmap('jet',22)
    caxis([0.0 22.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_toc_provinces,'-depsc2', ['./SWI_eps/Provinces_SouthPac.eps']);
    
    % plot a-value vs water-depth:
    % reshape a-values and water-depth to vextor
    
    a_value_Vec = reshape(Parameter_a_2Dtotal,[],1);
    water_depth_aligned_Vec = reshape(water_depth_aligned, [],1);
    
    figure;
    scatter(a_value_Vec,water_depth_aligned_Vec)
    
    
end

%% fill in the a-values for missing grid-cells with

