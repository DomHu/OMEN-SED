% plot locations where TOC matching at zBIO failed on top of TOC data
path(path,'/home/domhu/Documents/MATLAB/M_Map');


resolution = 1;     % 1: 1degree; 2: 2degree


switch resolution
    case 1
        
        load('./data/01_plot_BC_April_21/BC_1degree/Lee_toc_lr_weighted_SWI_PK_filled_a1to10_por085.mat');  % inversely calculated TOC at SWI [wt%]     -- need to translate to mol/cm^3 solid phase, i.e. *1e-2/12*bsd.rho_sed
        %                            load('./data/01_plot_BC_April_21/BC_1degree/Lee_toc_lr_weighted_SWI_a10to100_por085.mat');  % inversely calculated TOC at SWI [wt%]     -- need to translate to mol/cm^3 solid phase, i.e. *1e-2/12*bsd.rho_sed
        toc_in = TOC_SWI_BC;
        load('./data/01_plot_BC_April_21/BC_1degree/lat_lr.mat');
        load('./data/01_plot_BC_April_21/BC_1degree/long_lr.mat');
        long = long_lr;
        lat = lat_lr;
        
        ocean_cells = ~isnan(toc_in);
        Number_ocean_cells = sum(sum(ocean_cells));
        
        % get the coordinates for failed cells match at zbio
        Loc_zbiofail_lat = lat(debug.zbioMatching_xy_GMD(:,1));
        Loc_zbiofail_long = long(debug.zbioMatching_xy_GMD(:,2));
        
        % get the coordinates for failed cells - negative oxidation
        Loc_OXIDfail_lat = lat(debug.NEGATIVE_OXID_xy(:,1));
        Loc_OXIDfail_long = long(debug.NEGATIVE_OXID_xy(:,2));
        
    case 2
        
        load './data/BC_2degrees/toc_02.mat'                                       % TOC at SWI [wt%]     -- need to translate to mol/cm^3 solid phase, i.e. *1e-2/12*bsd.rho_sed
        toc_in = toc_02;
        load './data/BC_2degrees/lat_02.mat'
        load './data/BC_2degrees/long_02.mat'
        long = long_02;
        lat = lat_02;
        
        ocean_cells = ~isnan(toc_02);
        Number_ocean_cells = sum(sum(ocean_cells));
        
        % get the coordinates for failed cells
        Loc_zbiofail_lat = lat(debug.zbioMatching_xy(:,1));
        Loc_zbiofail_long = long(debug.zbioMatching_xy(:,2));
        
end

% using Robinson projection
fig_toc = figure;
m_proj('Robinson','longitudes',[-180 179.99], ...
    'latitudes',[-90 90]);
hold on;
% set(gca,'FontSize',30)
toc_levels = [0:0.1:15];
[C,h] = m_contourf(long, lat, toc_in, toc_levels);
set(h,'LineColor','none')
title('Seiter TOC (wt%) + failed TOC calculation')
m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
m_grid('fontsize',8);
m_plot(Loc_zbiofail_long, Loc_zbiofail_lat,'xr','MarkerSize',4);
m_plot(Loc_OXIDfail_long, Loc_OXIDfail_lat,'xg','MarkerSize',4);
formatSpec = '%.2f';
text1 = {['TOC-match fail: ' int2str(size(Loc_zbiofail_lat,1)), ' of ', int2str(Number_ocean_cells), ' cells']
        [ num2str(size(Loc_zbiofail_lat,1)/Number_ocean_cells*100, formatSpec), '%']};
text2 =  {       ['Cox < 0: ' int2str(size(Loc_OXIDfail_lat,1)), ' of ', int2str(Number_ocean_cells), ' cells']
        [ num2str(size(Loc_OXIDfail_lat,1)/Number_ocean_cells*100, formatSpec), '%']
        };
m_text(45, 60, text1,'FontSize',5, 'FontWeight', 'Bold', 'Color', 'r');
m_text(45, 48, text2,'FontSize',5, 'FontWeight', 'Bold', 'Color', 'g');
colorbar ('horizontal')
colormap(parula)
caxis([0.0 3.0])
xlabel('Longitude')
ylabel('Latitude')
print(fig_toc,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/999_OM_oxidation_FAILS_' str_date '.eps']);

