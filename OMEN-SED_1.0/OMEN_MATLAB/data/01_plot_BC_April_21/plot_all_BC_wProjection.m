
path(path,'/home/domhu/Documents/MATLAB/M_Map');

plot_BW_conditions = false;
plot_depth_sedrate = false;

plot_Lee_data = false;
plot_Lee_toc_calculatedSWI = true;  % plot the inversely calculated SWI values
plot_Seiter_data = false;
plot_a_value = false;
plot_no_projection = false;


if(plot_depth_sedrate)

      
% seafloor depth from NASA: https://neo.sci.gsfc.nasa.gov/view.php?datasetId=GEBCO_BATHY
load './data/GEBCO_BATHY_2002-01-01_rgb_1440x720.CSV'

water_depth = flipud(GEBCO_BATHY_2002_01_01_rgb_1440x720);
water_depth_nan=find(water_depth>=0);
water_depth(water_depth_nan)=NaN;

% % add row and column of NaNs - so it's same as toc data
% [n, m]=size(water_depth);
% row_nan = NaN(1,m);
% column_nan = NaN(1,n+1)';
% water_depth=[row_nan;water_depth];
% water_depth=[column_nan, water_depth];

% make lat / lon for 1/4 resolution
lat_quarter = [-89.875:0.25:89.875]';
long_quarter = [-179.875:0.25:179.875]';


% Robinson projection
fig_SFD = figure;
m_proj('Robinson','longitudes',[-180 179.99], ...
           'latitudes',[-90 90]);
hold on; 
% set(gca,'FontSize',30)
depth_levels=[-8000:200:0];
[C,h] = m_contourf(long_quarter, lat_quarter, water_depth,depth_levels);
set(h,'LineColor','none')
title('Water Depth (m)')
m_coast('linewidth',0.5); %('patch',[.7 .7 .7]);
m_grid('fontsize',8);
% hold on;
% shading interp; 
colorbar ('horizontal')
colormap(parula)
caxis([-5000.0 0.0])
xlabel('Longitude')
ylabel('Latitude')
print(fig_SFD,'-depsc2', ['SeafloorDepth_NASA_Robinson_HR.eps']);


% sedimentation rate, cm/yr (after Burwicz et al. (2011)) 
w1 = 0.117;
w2 = 0.006;
z1 = 200;
z2 = 4000;
c1 = 3;
c2 = 10;
sed_rate = w1./(1+(water_depth./z1).^c1) + w2./(1+(water_depth./z2).^c2);

fig_sedrate = figure;
m_proj('Robinson','longitudes',[-180 179.99], ...
           'latitudes',[-90 90]);
hold on; 
% set(gca,'FontSize',30)
    sed_levels = [0:0.005:0.13];
[C,h] = m_contourf(long_quarter, lat_quarter, sed_rate,sed_levels);
set(h,'LineColor','none')
    title('Sedimentation rate after Burwicz (cm/yr)')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
% hold on;
% shading interp; 
colorbar ('horizontal')
colormap(parula)
    caxis([0.0 0.12])
xlabel('Longitude')
ylabel('Latitude')
print(fig_sedrate,'-depsc2', ['Sed_Rate_Burwicz_Robinson_HR.eps']);


    
    % sedimentation rate holocene from Bradley
    
    load('./data/Bradley/long.dat');
    load('./data/Bradley/lat.dat');
    load('./data/Bradley/sed_holo.mat');
    
    % using Robinson projection
    fig_sedrate_Brad = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    sed_levels = [0:0.005:0.13];
    [C,h] = m_contourf(long, lat, sed_holo, sed_levels);    % 1/4 degree resolution
    set(h,'LineColor','none')
    title('Sedimentation rate as Bradley (cm/yr)')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 0.12])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_sedrate_Brad,'-depsc2', ['Sed_holocene_asBradley_Robinson.eps']);
    
    
end

if(plot_BW_conditions)
    
    % load data from WOA in 1 degree resolution
  	load('./data/WOA_2018_1degree/lat_WOA_01.mat');
  	load('./data/WOA_2018_1degree/long_WOA_01.mat');
  	load('./data/WOA_2018_1degree/Tmp_BW_2018_01.mat');
  	load('./data/WOA_2018_1degree/O2_BW_2018.mat');
  	load('./data/WOA_2018_1degree/NO3_BW_2018.mat');

  	% Robinson projection
    fig_tmp_BW = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    tmp_levels = [0:1:31];
    [C,h] = m_contourf(long_WOA_01, lat_WOA_01, Tmp_BW_2018_01, tmp_levels);
    set(h,'LineColor','none')
    title('Temperature BW (degree C)')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 30.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_tmp_BW,'-depsc2', ['Tmp_BW_Robinson.eps']);
    
    
    fig_O2_BW = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    O2_levels = [0:10:400];
    [C,h] = m_contourf(long_WOA_01, lat_WOA_01, O2_BW_2018, O2_levels);
    set(h,'LineColor','none')
    title('Oxygen BW (\muM)')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 350.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_O2_BW,'-depsc2', ['O2_BW_Robinson.eps']);
      
    fig_O2_BW = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    NO3_levels = [0:1:52];
    [C,h] = m_contourf(long_WOA_01, lat_WOA_01, NO3_BW_2018, NO3_levels);
    set(h,'LineColor','none')
    title('Nitrate BW (\muM)')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 50.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_O2_BW,'-depsc2', ['NO3_BW_Robinson.eps']);   
    
    
    
    % FeOOH calculated following Dale using sedimentation rate
    
    load('./data/FluxFe_total_Burw_BW_01_por85.mat');               % [mumol/(cm2 yr)]
    load('./data/FluxFe_total_Burw_Daleunits_BW_01_por85.mat');     % [mumol/(m2 d)]
        
    fig_O2_BW = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    feooh_levels = [0:25:1150];
    [C,h] = m_contourf(long_WOA_01, lat_WOA_01, FluxFe_total_Burw_Daleunits_BW_01_por85, feooh_levels);
    set(h,'LineColor','none')
    title('Total flux Fe(OH)3 BW (\mumol m-2 d-1)')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 800.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_O2_BW,'-depsc2', ['FeOOH_BW_Robinson.eps']);             
    
    
  

end

if(plot_Lee_toc_calculatedSWI)
       %% plot the inversely calculated SWI TOC wt%
 	load('./BC_1degree/lat_lr.mat');      
    load('./BC_1degree/long_lr.mat');      
  	load('./BC_1degree/Lee_toc_lr_weighted_SWI_a10to100_por085.mat');
    TOC_SWI_BC_10to100 = TOC_SWI_BC;
  	load('./BC_1degree/Lee_toc_lr_weighted_SWI_a1to100_por085.mat');      
    TOC_SWI_BC_1to100 = TOC_SWI_BC;
  	load('./BC_1degree/Lee_toc_lr_weighted.mat');      % original

  	% Robinson projection
    fig_toc_Lee = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    toc_levels = [-0.1:0.005:0.1];
    [C,h] = m_contourf(long_lr, lat_lr, TOC_SWI_BC_10to100 - Lee_toc_lr_weighted,toc_levels);
    set(h,'LineColor','none')
    title('TOC anomaly: SWI(calc) - Lee(mean) (wt%); a = 10-100')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([-0.1 0.1])
    xlabel('Longitude')
    ylabel('Latitude')
    
      	% Robinson projection
    fig_toc_Lee = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    toc_levels = [-0.1:0.005:0.1];
    [C,h] = m_contourf(long_lr, lat_lr, TOC_SWI_BC_1to100 - Lee_toc_lr_weighted,toc_levels);
    set(h,'LineColor','none')
    title('TOC anomaly: SWI(calc) - Lee(mean) (wt%); a = 1-100')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([-0.1 0.1])
    xlabel('Longitude')
    ylabel('Latitude')
    
      	% Robinson projection
    fig_toc_Lee = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    toc_levels = [0:0.1:5];
    [C,h] = m_contourf(long_lr, lat_lr, Lee_toc_lr_weighted,toc_levels);
    set(h,'LineColor','none')
    title('TOC original Lee(mean) (wt%)')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([0 3.0])
    xlabel('Longitude')
    ylabel('Latitude')    
        
        
      	% Robinson projection
    fig_toc_Lee = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    toc_levels = [0:0.1:5];
    [C,h] = m_contourf(long_lr, lat_lr, TOC_SWI_BC_1to100,toc_levels);
    set(h,'LineColor','none')
    title('TOC SWI(calc) (wt%); a = 1-100')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([0 3.0])
    xlabel('Longitude')
    ylabel('Latitude')
     
end

if(plot_Lee_data)
    %% load toc from Lee ea 2019
    load('./data/DatasetS1_Lee_ea2019.txt');
    % reshape into matrix
    Lee_TOC=DatasetS1_Lee_ea2019(:,3);
    Lee_TOC_matrix = reshape(Lee_TOC,4320,2160);
    Lee_long= DatasetS1_Lee_ea2019(1:4320,1);
    Lee_lat= DatasetS1_Lee_ea2019(1:4320:end,2);
    
    % Robinson projection
    fig_toc_Lee = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    toc_levels = [0:0.1:5];
    [C,h] = m_contourf(Lee_long, Lee_lat, Lee_TOC_matrix',toc_levels);
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
    print(fig_toc_Lee,'-depsc2', ['TOC_input_Lee_Robinson.eps']);
    
    % without projection
    if(plot_no_projection)
        fig_toc_Lee_noproj = figure;
        set(gca,'FontSize',30)
        pcolor(Lee_long, Lee_lat, Lee_TOC_matrix');
        title('Lee SWI-TOC (wt%)')
        hold on;
        shading interp;
        %contour(long, lat, toc,'LineColor','k')
        colorbar ('horizontal')
        colormap(parula)
        caxis([0.0 3.0])
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_toc_Lee_noproj,'-depsc2', ['TOC_input_Lee_noproj.eps']);
    end
    
end

if(plot_Seiter_data)
    
    %% plot TOC wt% as used in Bradley
    
    % 1/4 resolution
    load('./data/Bradley/toc.dat');
    load('./data/Bradley/long.dat');
    load('./data/Bradley/lat.dat');
    load('./data/Bradley/sed_holo.mat');
    
    TOC_nan=find(toc<0);
    toc(TOC_nan)=NaN;
    % save('toc_NaN.mat' , 'toc')
    
    % % 1 degree
    % load '../BC_1degrees/toc_01.mat'                                       % TOC at SWI [wt%]     -- need to translate to mol/cm^3 solid phase, i.e. *1e-2/12*bsd.rho_sed
    % toc_in = toc_01;
    % % long & lat from WOA in 1 degree resolution
    % load('../lat_WOA_01.mat')
    % load('../long_WOA_01.mat')
    
    % using Robinson projection
    fig_toc = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    % m_contourf(long_WOA_01, lat_WOA_01, toc_in);    % 1 degree resolution
    toc_levels = [0:0.1:15];
    [C,h] = m_contourf(long, lat, toc, toc_levels);    % 1/4 degree resolution
    set(h,'LineColor','none')
    title('Seiter TOC as in Bradley (wt%)')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 3.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_toc,'-depsc2', ['TOC_input_Seiter_Robinson.eps']);
    

    
    %% Pika TOC
    load('./data/Pika/Parameter_log10_a_with_BC.mat');
    % lon = T_1D_Grids(:,1);
    % lat = T_1D_Grids(:,2);
    
    % convert into array
    Array_BC = table2array(T_1D_Grids);
    
    TOC_Seiter = Array_BC(:,3);
    
    % reshape TOC into matrix
    TOC_Seiter_matrix = reshape(TOC_Seiter,180,360);
    Seiter_long= Array_BC(1:180:end,1);
    Seiter_lat= Array_BC(1:180,2);
    
    fig_toc_Seiter = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    toc_levels = [0:0.1:11];
    [C,h] = m_contourf(Seiter_long, Seiter_lat, TOC_Seiter_matrix,  toc_levels);
    set(h,'LineColor','none')
    title('Seiter SWI-TOC PK table (wt%)')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 3.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_toc_Seiter,'-depsc2', ['TOC_input_Seiter_Robinson_PKtable.eps']);
    
    
    
    % without projection
    if(plot_no_projection)
        fig_SWI_TOC = figure;
        set(gca,'FontSize',30)
        pcolor(long, lat, toc);
        title('Seiter TOC as in Bradley (wt%)')
        hold on;
        shading interp;
        colorbar ('horizontal')
        colormap(parula)
        caxis([0.0 3.0])
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_SWI_TOC,'-depsc2', ['TOC_input_Seiter_noproj.eps']);
    end
    
   
end


if(plot_a_value)
    
    %% Parameter a in 2D  -- this is in 1 degree resolution
    load('./data/Pika/Parameter_log10_a_2D.mat')
    % long & lat from WOA in 1 degree resolution
    load('./data/WOA_2018_1degree/lat_WOA_01.mat')
    load('./data/WOA_2018_1degree/long_WOA_01.mat')
    
    Parameter_a_2D = flipud(Parameter_a_2D);
    Parameter_a_2Dtotal = 10.^Parameter_a_2D;
    
    fig01 = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    a_levels = [-0.1:0.1:2.0];
    [C,h] = m_contourf(long_WOA_01, lat_WOA_01, Parameter_a_2D);
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
    
    
    
end

