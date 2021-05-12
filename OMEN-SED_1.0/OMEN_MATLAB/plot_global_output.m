% function plot_global_output(high_res)
%% plots the results of the global analysis for OM degradation pathways
% high_res:     true    -   1/4 degree resolution
%               false   -   1 degree resolution

print_SWI_TOC_BC = false;
calc_TOC_burial = false;
calc_Regional_OM_oxidation = false;
calc_Total_OM_oxidation = false;    % Global total in mol/yr -- taking size of grid-cells into account!
plot_oxidation_rates = false;
plot_penetration = false;
plot_SWI_fluxes = false;
calc_plot_hyp_remin_rates = true;

plot_BC = false;

% NOT NEEDED:
plot_Fe2_efflux = false;
plot_oxidation_rates_fluxes = false;  %% not needed - is the same as calculate rates with calcReac (plot_oxidation_rates)

zinf = 100;

% high_res = false;
calc_res = 2;           % 1: 1/4°;  2: 1°;  3: 2°

switch calc_res
    case 1
        load data/long.dat
        load data/lat.dat
        load './data/water_depth_updated.mat'                           % Seafloor depth [m]
        water_depth_updated = -water_depth_updated;
        
    case 2
        load data/BC_1degrees/long_WOA_01.mat
        load data/BC_1degrees/lat_WOA_01.mat
        long = long_WOA_01;
        lat = lat_WOA_01;
        load('./data/01_plot_BC_April_21/BC_1degree/water_depth_aligned.mat');                             % Seafloor depth [m]
        water_depth_updated = -water_depth_aligned;

        load('./data/01_plot_BC_April_21/BC_1degree/Lee_toc_lr_weighted_SWI_PK_filled_a1to10_por085.mat');  % inversely calculated TOC at SWI [wt%]     -- need to translate to mol/cm^3 solid phase, i.e. *1e-2/12*bsd.rho_sed
        toc_in = TOC_SWI_BC;
        
        
    case 3
        load data/BC_2degrees/long_02.mat
        load data/BC_2degrees/lat_02.mat
        long = long_02;
        lat = lat_02;
        load './data/BC_2degrees/water_depth_02_updated.mat'                           % Seafloor depth [m]
        water_depth_updated = water_depth_02_updated;
        
end


str_date = [datestr(date,7), datestr(date,5), datestr(date,11)];
lightGrey2   = [0.85 0.85 0.85];

if(print_SWI_TOC_BC)
%  PRINT THE TOC SWI concentration
    fig_SWI_fluxO2 = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [0:0.01:5.0];
        limits = [0 3.0];
        [C,h] = m_contourf(long, lat,  toc_in, levels);
        set(h,'LineColor','none')
        title('TOC (wt%) -- Lee et al. (2019)')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_SWI_fluxO2,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/SWI_TOC_inverse_' str_date '.eps']);    
end

if(calc_TOC_burial)
    % in mol cm-2 yr-1
    %     load output/0302_50cm_100G_Arndt_rescales_1-10_1degree_TOC_burial_flux_030221.mat
    
    % define depth regimes after Dunne et al. 2005:
    near_shore = 50;
    shelf = 200;
    slope = 2000;
    
    %% Global burial values
    TOC_burial_flux_m2 = -TOC_burial_flux.*100^2;    % convert from mol cm-2 yr-1  to mol m-2 yr-1
    % Global burial rate in mol yr-1
    Total_OM_burial_mol_Global = nansum(nansum(TOC_burial_flux_m2.*dxdy));
    % Global burial rate in Pg yr-1 (from mmol yr-1)
    Total_OM_burial_Pg.Global = 12*(Total_OM_burial_mol_Global)*10^(-15);
    
    %% Depth regime burial values in Pg C yr-1
    
    Total_OM_burial_Pg.NearShore = nansum(nansum(TOC_burial_flux_m2(water_depth_updated<near_shore).*dxdy(water_depth_updated<near_shore)))*12*10^(-15);
    Total_OM_burial_Pg.Shelf = nansum(nansum(TOC_burial_flux_m2(water_depth_updated>=near_shore & water_depth_updated<shelf).*dxdy(water_depth_updated>=near_shore & water_depth_updated<shelf)))*12*10^(-15);
    Total_OM_burial_Pg.Slope = nansum(nansum(TOC_burial_flux_m2(water_depth_updated>=shelf & water_depth_updated<slope).*dxdy(water_depth_updated>=shelf & water_depth_updated<slope)))*12*10^(-15);
    Total_OM_burial_Pg.Plain = nansum(nansum(TOC_burial_flux_m2(water_depth_updated>=slope).*dxdy(water_depth_updated>=slope)))*12*10^(-15);
    
end


if(calc_Regional_OM_oxidation)
    % define depth regimes after Dunne et al. 2005:
    near_shore = 50;
    shelf = 200;
    slope = 2000;
    
    %     % m^2 per grid-cell
    %     load output/dxdy_0302_50cm_100G_Arndt_rescales_1-10_2degree_030221.mat
    %
    %     % oxidation rates in mmol m-2 yr-1
    %     load output/Cox_rate_out_0302_50cm_100G_Arndt_rescales_1-10_2degree_030221.mat
    Cox_aerobic = Cox_rate_out.Cox_aerobic;
    Cox_denitr=Cox_rate_out.Cox_denitr;
    Cox_FeIII= Cox_rate_out.Cox_FeIII;
    Cox_sulfate=Cox_rate_out.Cox_sulfate;
    Cox_total=Cox_rate_out.Cox_total;
    
    % calculate global oxid.rates in mmol yr-1
    Total_OM_oxid_mmol_Global = nansum(nansum(Cox_total.*dxdy));
    % oxidation rate in Pg yr-1 (from mmol yr-1)
    Total_OM_oxid_Pg.Global = 12*(Total_OM_oxid_mmol_Global/1000)*10^(-15);
    
    %% Total depth regime oxidation rates in mmol yr-1
    Total_OM_oxid_mmol_NearShore = nansum(nansum(Cox_total(water_depth_updated<near_shore).*dxdy(water_depth_updated<near_shore)));
    Total_OM_oxid_mmol_Shelf = nansum(nansum(Cox_total(water_depth_updated>=near_shore & water_depth_updated<shelf).*dxdy(water_depth_updated>=near_shore & water_depth_updated<shelf)));
    Total_OM_oxid_mmol_Slope = nansum(nansum(Cox_total(water_depth_updated>=shelf & water_depth_updated<slope).*dxdy(water_depth_updated>=shelf & water_depth_updated<slope)));
    Total_OM_oxid_mmol_Plain = nansum(nansum(Cox_total(water_depth_updated>=slope).*dxdy(water_depth_updated>=slope)));
    
    % Total depth regime oxidation rates in Pg yr-1 (from mmol yr-1)
    Total_OM_oxid_Pg.NearShore = 12*(Total_OM_oxid_mmol_NearShore/1000)*10^(-15);
    Total_OM_oxid_Pg.Shelf = 12*(Total_OM_oxid_mmol_Shelf/1000)*10^(-15);
    Total_OM_oxid_Pg.Slope = 12*(Total_OM_oxid_mmol_Slope/1000)*10^(-15);
    Total_OM_oxid_Pg.Plain = 12*(Total_OM_oxid_mmol_Plain/1000)*10^(-15);
    
    
    %% TEA specific depth regime oxidation rates in Pg C yr-1
    % Near shore
    NearShore_OM_oxid_Pg.aerobic = nansum(nansum(Cox_aerobic(water_depth_updated<near_shore).*dxdy(water_depth_updated<near_shore)))*12/1000*10^(-15);
    NearShore_OM_oxid_Pg.denitr = nansum(nansum(Cox_denitr(water_depth_updated<near_shore).*dxdy(water_depth_updated<near_shore)))*12/1000*10^(-15);
    NearShore_OM_oxid_Pg.FeIII = nansum(nansum(Cox_FeIII(water_depth_updated<near_shore).*dxdy(water_depth_updated<near_shore)))*12/1000*10^(-15);
    NearShore_OM_oxid_Pg.sulfate = nansum(nansum(Cox_sulfate(water_depth_updated<near_shore).*dxdy(water_depth_updated<near_shore)))*12/1000*10^(-15);
    
    Frac_NearShore.aerobic = NearShore_OM_oxid_Pg.aerobic/Total_OM_oxid_Pg.NearShore*100;
    Frac_NearShore.denitr = NearShore_OM_oxid_Pg.denitr/Total_OM_oxid_Pg.NearShore*100;
    Frac_NearShore.FeIII = NearShore_OM_oxid_Pg.FeIII/Total_OM_oxid_Pg.NearShore*100;
    Frac_NearShore.sulfate = NearShore_OM_oxid_Pg.sulfate/Total_OM_oxid_Pg.NearShore*100;
    
    if(abs(NearShore_OM_oxid_Pg.aerobic+NearShore_OM_oxid_Pg.denitr+NearShore_OM_oxid_Pg.FeIII+NearShore_OM_oxid_Pg.sulfate - Total_OM_oxid_Pg.NearShore)>1e-6)
        error('NearShore ~= 1.0');
    end
    
    % Shelf
    Shelf_OM_oxid_Pg.aerobic = nansum(nansum(Cox_aerobic(water_depth_updated>=near_shore & water_depth_updated<shelf).*dxdy(water_depth_updated>=near_shore & water_depth_updated<shelf)))*12/1000*10^(-15);
    Shelf_OM_oxid_Pg.denitr = nansum(nansum(Cox_denitr(water_depth_updated>=near_shore & water_depth_updated<shelf).*dxdy(water_depth_updated>=near_shore & water_depth_updated<shelf)))*12/1000*10^(-15);
    Shelf_OM_oxid_Pg.FeIII = nansum(nansum(Cox_FeIII(water_depth_updated>=near_shore & water_depth_updated<shelf).*dxdy(water_depth_updated>=near_shore & water_depth_updated<shelf)))*12/1000*10^(-15);
    Shelf_OM_oxid_Pg.sulfate = nansum(nansum(Cox_sulfate(water_depth_updated>=near_shore & water_depth_updated<shelf).*dxdy(water_depth_updated>=near_shore & water_depth_updated<shelf)))*12/1000*10^(-15);
    
    Frac_Shelf.aerobic = Shelf_OM_oxid_Pg.aerobic/Total_OM_oxid_Pg.Shelf*100;
    Frac_Shelf.denitr = Shelf_OM_oxid_Pg.denitr/Total_OM_oxid_Pg.Shelf*100;
    Frac_Shelf.FeIII = Shelf_OM_oxid_Pg.FeIII/Total_OM_oxid_Pg.Shelf*100;
    Frac_Shelf.sulfate = Shelf_OM_oxid_Pg.sulfate/Total_OM_oxid_Pg.Shelf*100;
    
    if(abs(Shelf_OM_oxid_Pg.aerobic+Shelf_OM_oxid_Pg.denitr+Shelf_OM_oxid_Pg.FeIII+Shelf_OM_oxid_Pg.sulfate - Total_OM_oxid_Pg.Shelf)>1e-6)
        error('Shelf ~= 1.0');
    end
    
    % Slope
    Slope_OM_oxid_Pg.aerobic = nansum(nansum(Cox_aerobic(water_depth_updated>=shelf & water_depth_updated<slope).*dxdy(water_depth_updated>=shelf & water_depth_updated<slope)))*12/1000*10^(-15);
    Slope_OM_oxid_Pg.denitr = nansum(nansum(Cox_denitr(water_depth_updated>=shelf & water_depth_updated<slope).*dxdy(water_depth_updated>=shelf & water_depth_updated<slope)))*12/1000*10^(-15);
    Slope_OM_oxid_Pg.FeIII = nansum(nansum(Cox_FeIII(water_depth_updated>=shelf & water_depth_updated<slope).*dxdy(water_depth_updated>=shelf & water_depth_updated<slope)))*12/1000*10^(-15);
    Slope_OM_oxid_Pg.sulfate = nansum(nansum(Cox_sulfate(water_depth_updated>=shelf & water_depth_updated<slope).*dxdy(water_depth_updated>=shelf & water_depth_updated<slope)))*12/1000*10^(-15);
    
    Frac_Slope.aerobic = Slope_OM_oxid_Pg.aerobic/Total_OM_oxid_Pg.Slope*100;
    Frac_Slope.denitr = Slope_OM_oxid_Pg.denitr/Total_OM_oxid_Pg.Slope*100;
    Frac_Slope.FeIII = Slope_OM_oxid_Pg.FeIII/Total_OM_oxid_Pg.Slope*100;
    Frac_Slope.sulfate = Slope_OM_oxid_Pg.sulfate/Total_OM_oxid_Pg.Slope*100;
    
    if(abs(Slope_OM_oxid_Pg.aerobic+Slope_OM_oxid_Pg.denitr+Slope_OM_oxid_Pg.FeIII+Slope_OM_oxid_Pg.sulfate - Total_OM_oxid_Pg.Slope)>1e-6)
        error('Slope ~= 1.0');
    end
    
    % Plain
    Plain_OM_oxid_Pg.aerobic = nansum(nansum(Cox_aerobic(water_depth_updated>=slope).*dxdy(water_depth_updated>=slope)))*12/1000*10^(-15);
    Plain_OM_oxid_Pg.denitr = nansum(nansum(Cox_denitr(water_depth_updated>=slope).*dxdy(water_depth_updated>=slope)))*12/1000*10^(-15);
    Plain_OM_oxid_Pg.FeIII = nansum(nansum(Cox_FeIII(water_depth_updated>=slope).*dxdy(water_depth_updated>=slope)))*12/1000*10^(-15);
    Plain_OM_oxid_Pg.sulfate = nansum(nansum(Cox_sulfate(water_depth_updated>=slope).*dxdy(water_depth_updated>=slope)))*12/1000*10^(-15);
    
    Frac_Plain.aerobic = Plain_OM_oxid_Pg.aerobic/Total_OM_oxid_Pg.Plain*100;
    Frac_Plain.denitr = Plain_OM_oxid_Pg.denitr/Total_OM_oxid_Pg.Plain*100;
    Frac_Plain.FeIII = Plain_OM_oxid_Pg.FeIII/Total_OM_oxid_Pg.Plain*100;
    Frac_Plain.sulfate = Plain_OM_oxid_Pg.sulfate/Total_OM_oxid_Pg.Plain*100;
    
    if(abs(Plain_OM_oxid_Pg.aerobic+Plain_OM_oxid_Pg.denitr+Plain_OM_oxid_Pg.FeIII+Plain_OM_oxid_Pg.sulfate - Total_OM_oxid_Pg.Plain)>1e-6)
        error('Plain ~= 1.0');
    end
    
    %% calculate area of depth regimes
    Area.NearShore = nansum(nansum(dxdy(water_depth_updated<near_shore)));
    Area.Shelf = nansum(nansum(dxdy(water_depth_updated>=near_shore & water_depth_updated<shelf)));
    Area.Slope = nansum(nansum(dxdy(water_depth_updated>=shelf & water_depth_updated<slope)));
    Area.Plain = nansum(nansum(dxdy(water_depth_updated>=slope)));
    
end

if(calc_Total_OM_oxidation)
    
    % load output/0502_100cm_100G_Arndt_rescales_1-10_1degree_gamma0.95_dxdy_050221.mat
    %
    %      % oxidation rates in mmol m-2 yr-1
    %    load output/Cox_rate_out_2m_1xTOC_gammaFe2_par_Fe3por0.85_SMIsinkDale_gamma0.95_101220.mat
    Cox_aerobic = Cox_rate_out.Cox_aerobic;
    Cox_denitr=Cox_rate_out.Cox_denitr;
    Cox_FeIII= Cox_rate_out.Cox_FeIII;
    Cox_sulfate=Cox_rate_out.Cox_sulfate;
    Cox_total=Cox_rate_out.Cox_total;
    
    % calculate global oxid.rates in mmol yr-1
    Global_OM_oxid_mmol.Total = nansum(nansum(Cox_total.*dxdy));
    Global_OM_oxid_Tmolyr = (Global_OM_oxid_mmol.Total/1000)/10^(12);
    Global_OM_oxid_mmol.aerobic = nansum(nansum(Cox_aerobic.*dxdy));
    Global_OM_oxid_mmol.denitr = nansum(nansum(Cox_denitr.*dxdy));
    Global_OM_oxid_mmol.FeIII = nansum(nansum(Cox_FeIII.*dxdy));
    Global_OM_oxid_mmol.sulfate = nansum(nansum(Cox_sulfate.*dxdy));
    
    Global_OM_oxid_Frac.aerobic = Global_OM_oxid_mmol.aerobic/Global_OM_oxid_mmol.Total*100;
    Global_OM_oxid_Frac.denitr = Global_OM_oxid_mmol.denitr/Global_OM_oxid_mmol.Total*100;
    Global_OM_oxid_Frac.FeIII = Global_OM_oxid_mmol.FeIII/Global_OM_oxid_mmol.Total*100;
    Global_OM_oxid_Frac.sulfate = Global_OM_oxid_mmol.sulfate/Global_OM_oxid_mmol.Total*100;
    
    % oxidation rate in Pg yr-1 (from mmol yr-1)
    Global_OM_oxid_Total_Pg = 12*(Global_OM_oxid_mmol.Total/1000)*10^(-15);
    
    %     fig_dxdy= figure;
    %     %set(gca,'FontSize',30)
    %     pcolor(long, lat, dxdy);
    %     hold on;
    %     % set(gca,'Color',lightGrey2)
    %     % set(gcf,'inverthardcopy','off');
    %     title('Size of grid cell (m^2)')
    %     shading interp;
    %     %contour(long, lat, toc,'LineColor','k')
    %     colorbar ('horizontal')
    %     colormap(parula)
    % %    caxis([0 12])
    %     xlabel('Longitude')
    %     ylabel('Latitude')
    %     print(fig_dxdy,'-depsc2', ['./plots/' exp_name '_dxdy_grid-size_' str_date '.eps']);
    %
end

%% plot oxidation rates total and fraction per pathway; calculated as reaction rates in layers -- calcReac (saved in mmol m-2 day-1)
if(plot_oxidation_rates)
    % oxidation rates in mmol m-2 yr-1
    % load output/0502_100cm_100G_Arndt_rescales_1-10_1degree_gamma0.95_Cox_rate_out_050221.mat
    Cox_aerobic = Cox_rate_out.Cox_aerobic;
    Cox_denitr=Cox_rate_out.Cox_denitr;
    Cox_FeIII= Cox_rate_out.Cox_FeIII;
    Cox_sulfate=Cox_rate_out.Cox_sulfate;
    Cox_total=Cox_rate_out.Cox_total;
    Total_rate_sum=Cox_rate_out.Total_rate_sum;
    
    Cox_total_mumol_cmyr = Cox_total/1000*100^(-2)*10^6;    % from mmol m-2 yr-1 --> mumol cm-2 yr-1
    
    fig_OMoxid_total = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:1:250];
    limits = [0 200];
    [C,h] = m_contourf(long, lat, Cox_total_mumol_cmyr, levels);
    set(h,'LineColor','none')
    title('Total OM oxidation (\mumol cm-2 yr-1)')
    %        m_coast('linewidth',1,'color','k'); %
    m_coast('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_OMoxid_total,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/OM_oxidation_Total_' str_date '.eps']);
    
    
    % rates for diff TEA in mumol cm-2 yr-1
    fig_OM_ox_aerobic= figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:1:50];
    limits = [0 20];
    [C,h] = m_contourf(long, lat,  (Cox_aerobic/1000)*100^(-2)*10^6, levels);
    set(h,'LineColor','none')
    title('Aerobic degradation (\mumol cm-2 yr-1)')
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_OM_ox_aerobic,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/OM_oxidation_aerobic_' str_date '.eps']);
    
    
    fig_OM_ox_denit= figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:1:50];
    limits = [0 20];
    [C,h] = m_contourf(long, lat,  (Cox_denitr/1000)*100^(-2)*10^6, levels);
    set(h,'LineColor','none')
    title('Denitrification  (\mumol cm-2 yr-1)')
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_OM_ox_denit,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/OM_oxidation_denit_' str_date '.eps']);
    
    
    fig_OM_ox_Fe_red= figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:1:10];
    limits = [0 5];
    [C,h] = m_contourf(long, lat,  (Cox_FeIII/1000)*100^(-2)*10^6, levels);
    set(h,'LineColor','none')
    title('Fe-reduction  (\mumol cm-2 yr-1)')
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_OM_ox_Fe_red,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/OM_oxidation_FeIII_red_' str_date '.eps']);
    
    
    fig_OM_ox_sulfate_red= figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:1:150];
    limits = [0 100];
    [C,h] = m_contourf(long, lat,  (Cox_sulfate/1000)*100^(-2)*10^6, levels);
    set(h,'LineColor','none')
    title('Sulfate-reduction (in \mumol cm-2 yr-1)')
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_OM_ox_sulfate_red,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/OM_oxidation_Sulfate_red_' str_date '.eps']);
    
    if(true)
        %% Fraction of total
        fig_OM_ox_aerobic_frac= figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [0:0.5:100];
        limits = [0 100];
        [C,h] = m_contourf(long, lat,  Cox_aerobic./Cox_total*100, levels);
        set(h,'LineColor','none')
        title('Fraction of aerobic degradation (in %)')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_OM_ox_aerobic_frac,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/OM_oxidation_aerobic_red_frac_' str_date '.eps']);
        
        
        fig_OM_ox_denit_frac= figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [0:0.5:100];
        limits = [0 50];
        [C,h] = m_contourf(long, lat,  Cox_denitr./Cox_total*100, levels);
        set(h,'LineColor','none')
        title('Fraction of denitrification (in %)')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_OM_ox_denit_frac,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/OM_oxidation_denit_frac_' str_date '.eps']);
        
        
        fig_OM_ox_Fe_frac= figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [0:0.5:100];
        limits = [0 10];
        [C,h] = m_contourf(long, lat,  Cox_FeIII./Cox_total*100, levels);
        set(h,'LineColor','none')
        title('Fraction of Fe-reduction (in %)')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_OM_ox_Fe_frac,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/OM_oxidation_FeIII_red_frac_' str_date '.eps']);
        
        fig_OM_ox_sulfate_frac= figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [0:0.5:100];
        limits = [0 100];
        [C,h] = m_contourf(long, lat,  Cox_sulfate./Cox_total*100, levels);
        set(h,'LineColor','none')
        title('Fraction of sulfate-reduction (in %)')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_OM_ox_sulfate_frac,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/OM_oxidation_sulfate_red_frac_' str_date '.eps']);
        
        
    end
end




%% plot TEA penetration depth
if(plot_penetration)
    % load output/0502_100cm_100G_Arndt_rescales_1-10_1degree_gamma0.95_Penetration_out_050221.mat
    zox = Penetration_out.zox;
    zno3 = Penetration_out.zno3;
    zfeIII = Penetration_out.zfeIII;
    zso4 = Penetration_out.zso4;
    
    % if no penetration depth plot as NaN
    % zno3_inf = zno3;
    % zno3_inf(:,:) = NaN;
    zno3(zox==zinf) = NaN;
    zfeIII(zox==zinf) = NaN;
    zfeIII(zno3==zinf) = NaN;
    zso4(zox==zinf) = NaN;
    zso4(zno3==zinf) = NaN;
    zso4(zfeIII==zinf)  = NaN;
    
    
    fig_zox = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [0:0.5:100];
        limits = [0 100];
        [C,h] = m_contourf(long, lat,  zox, levels);
        set(h,'LineColor','none')
    title('O_2 penetration depth (cm)')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_zox,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/zox_' str_date '.eps']);

        
    fig_zNO3 = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [0:0.5:100];
        limits = [0 10];
        [C,h] = m_contourf(long, lat,  zno3, levels);
        set(h,'LineColor','none')
    title('NO_3 penetration depth (cm)')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_zNO3,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/zNO3_' str_date '.eps']);

   
        fig_zFeIII = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [0:0.5:100];
        limits = [0 10];
        [C,h] = m_contourf(long, lat,  zfeIII, levels);
        set(h,'LineColor','none')
        title('Fe-oxide penetration depth (cm)')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_zFeIII,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/zFeIII_' str_date '.eps']);

    
        fig_zSO4 = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [0:0.5:100];
        limits = [0 100];
        [C,h] = m_contourf(long, lat,  zso4, levels);
        set(h,'LineColor','none')
        title('SO_4 penetration depth (cm)')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_zSO4,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/zSO4_' str_date '.eps']);

end

%% plot calculated SWI-fluxes in mumol cm-2 yr-1
% positive values are fluxes into the sediments!
if(plot_SWI_fluxes)
    
    % load output/0502_100cm_100G_Arndt_rescales_1-10_1degree_gamma0.95_SWI_fluxes_050221.mat
    % load output/0502_100cm_100G_Arndt_rescales_1-10_1degree_gamma0.95_Flux_Fe2_Dale_units_050221.mat
    
 if(false)   
    fig_SWI_fluxO2 = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [0:0.5:200];
        limits = [0 100];
        [C,h] = m_contourf(long, lat,  -SWI_fluxes.flxswiO2*10^6, levels);
        set(h,'LineColor','none')
        title('SWI-flux O_2 (\mumol cm-2 yr-1)')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_SWI_fluxO2,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/SWI_flux_O2_' str_date '.eps']);    
    

    
    fig_SWI_fluxNO3 = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [-20:0.5:20];
        limits = [-10 10];
        [C,h] = m_contourf(long, lat,  -SWI_fluxes.flxswiNO3*10^6, levels);
        set(h,'LineColor','none')
        title('SWI-flux NO_3 (\mumol cm-2 yr-1)')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_SWI_fluxNO3,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/SWI_flux_NO3_' str_date '.eps']);    
    
    
    fig_SWI_fluxNH4 = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [-10:0.5:0];
        limits = [-2 0];
        [C,h] = m_contourf(long, lat,  -SWI_fluxes.flxswiNH4*10^6, levels);
        set(h,'LineColor','none')
        title('SWI-flux NH_4 (\mumol cm-2 yr-1)')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        colormap(flipud(parula))
       caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_SWI_fluxNH4,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/SWI_flux_NH4_' str_date '.eps']);    

    % units as Dale ea map Fig. 6a+c log(mmol m-2 yr-1)
    Flux_Fe2_logmmol = log10(Flux_Fe2_Dale_units*365*10^(-3));
    Flux_Fe2_logmmol(Flux_Fe2_logmmol==-Inf)=-10;
        Flux_Fe2_logmmol(Flux_Fe2_logmmol<-10)=-10;     % crazy low values set to -10 so they get plotted

    fig_Fe2_efflux = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [-10.5:0.05:2.0];
        limits = [-4 1.0];
        [C,h] = m_contourf(long, lat,  Flux_Fe2_logmmol, levels);
        set(h,'LineColor','none')
        title('Fe^{2+} SWI-flux (log(mmol m-2 yr-1))')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
       caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_Fe2_efflux,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/SWI_flux_Fe2_logmmolm2yr_' str_date '.eps']);    
        
   
    
    Flux_Fe2_logmumolcm = log10(Flux_Fe2_Dale_units*365*100^(-2));  % Dale units are umol m-2 d-1 as in Dale ea. 2015, Fig. 3 -- so translate to umol cm-2 yr-1
    Flux_Fe2_logmumolcm(Flux_Fe2_logmumolcm==-Inf)=-10;
   Flux_Fe2_logmumolcm(Flux_Fe2_logmumolcm<-10)=-10;     % crazy low values set to -10 so they get plotted

    fig_Fe2_efflux_myunits = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [-10.5:0.05:2.0];
        limits = [-4 1.0];
        [C,h] = m_contourf(long, lat,  Flux_Fe2_logmumolcm, levels);
        set(h,'LineColor','none')
    title('Fe^{2+} SWI-flux (log(\mumol cm-2 yr-1))')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
       caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_Fe2_efflux_myunits,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/SWI_flux_Fe2_logmumolcm2yr_' str_date '.eps']);    
    
end
        
    fig_SWI_fluxSO4 = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [0:0.5:100];
        limits = [0 10];
        [C,h] = m_contourf(long, lat,  -SWI_fluxes.flxswiSO4*10^6, levels);
        set(h,'LineColor','none')
    title('SWI-flux SO_4 (\mumol cm-2 yr-1)')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_SWI_fluxSO4,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/SWI_flux_SO4_' str_date '.eps']);    
        
      

    fig_SWI_fluxH2S = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [-20:0.5:0.0];
        limits = [-10 0];
        [C,h] = m_contourf(long, lat,  -SWI_fluxes.flxswiH2S*10^6, levels);
        set(h,'LineColor','none')
        title('SWI-flux H_2S (\mumol cm-2 yr-1)')
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        colormap(flipud(parula))
       caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_SWI_fluxH2S,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/SWI_flux_H2S_' str_date '.eps']);    

 

    
end


if(calc_plot_hyp_remin_rates)
    %% plot global hypsometry of remineralization rates (compare Fig. 5a-e, Thullner ea. 2009) and penetration depths
        % Gamma = 0.95
%     load output/0505_Lee_PK_filled_a1to10_GMDandSAIntConst_1degree_Inverse_Cox_rate_out_050521.mat
%     load output/0505_Lee_PK_filled_a1to10_GMDandSAIntConst_1degree_Inverse_dxdy_050521.mat
%     load output/0505_Lee_PK_filled_a1to10_GMDandSAIntConst_1degree_Inverse_SWI_fluxes_050521.mat
%     load output/0505_Lee_PK_filled_a1to10_GMDandSAIntConst_1degree_Inverse_Flux_Fe2_Dale_units_050521.mat
%     load output/0505_Lee_PK_filled_a1to10_GMDandSAIntConst_1degree_Inverse_Penetration_out_050521.mat
    
    % Gamma = 0.5
    load output/0605_Lee_PK_filled_a1to10_GMDandSAIntConst_1degree_Inverse_GAMMA05_Cox_rate_out_060521.mat
    load output/0605_Lee_PK_filled_a1to10_GMDandSAIntConst_1degree_Inverse_GAMMA05_dxdy_060521.mat
    load output/0605_Lee_PK_filled_a1to10_GMDandSAIntConst_1degree_Inverse_GAMMA05_SWI_fluxes_060521.mat
    load output/0605_Lee_PK_filled_a1to10_GMDandSAIntConst_1degree_Inverse_GAMMA05_Flux_Fe2_Dale_units_060521.mat
    load output/0605_Lee_PK_filled_a1to10_GMDandSAIntConst_1degree_Inverse_GAMMA05_Penetration_out_060521.mat
    
    hypsometry_depths = [0 50 150 350 750 1500 2750 4250 10000];
    %% calculate area of depth regimes
    for i = 2:length(hypsometry_depths)
        hypsometry_Area(i-1) = nansum(nansum(dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))));
        hypsometry_meanDepth(i-1) = nansum(nansum(dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*water_depth_updated(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/hypsometry_Area(i-1);
    end
    
    
    
    %     % oxidation rates in mmol m-2 yr-1
    %     load output/Cox_rate_out_0302_50cm_100G_Arndt_rescales_1-10_2degree_030221.mat
    Cox_aerobic = Cox_rate_out.Cox_aerobic;
    Cox_denitr=Cox_rate_out.Cox_denitr;
    Cox_FeIII= Cox_rate_out.Cox_FeIII;
    Cox_sulfate=Cox_rate_out.Cox_sulfate;
    Cox_total=Cox_rate_out.Cox_total;
    
    % SWI-fluxes all in mumol cm-2 yr-1  -- negative represents "out of sediment flux"
    SWIflux_O2 = -SWI_fluxes.flxswiO2*10^6;
    SWIflux_NO3 = -SWI_fluxes.flxswiNO3*10^6;
    SWIflux_NH4 = -SWI_fluxes.flxswiNH4*10^6;
    SWIflux_FeII = -Flux_Fe2_Dale_units*365*100^(-2);    	% Dale units are umol m-2 d-1 as in Dale ea. 2015, Fig. 3 -- so translate to umol cm-2 yr-1
    SWIflux_SO4 = -SWI_fluxes.flxswiSO4*10^6;
    SWIflux_H2S = -SWI_fluxes.flxswiH2S*10^6;
    
    % Penetration depths in cm
    zox = Penetration_out.zox;
    zno3 = Penetration_out.zno3;
    zfeIII = Penetration_out.zfeIII;
    zso4 = Penetration_out.zso4;
    %     % if no penetration depth plot as NaN
    %     zno3(zox==zinf) = NaN;
    %     zfeIII(zox==zinf) = NaN;
    %     zfeIII(zno3==zinf) = NaN;
    %     zso4(zox==zinf) = NaN;
    %     zso4(zno3==zinf) = NaN;
    %     zso4(zfeIII==zinf)  = NaN;
    
    for i = 2:length(hypsometry_depths)
        % Oxidation rates for the hysometry in mol m-2 yr-1
        Hypsometry_oxid_rates_aerobic_molperm2(i-1) = nansum(nansum(Cox_aerobic(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(1000*hypsometry_Area(i-1));
        Hypsometry_oxid_rates_denitr_molperm2(i-1) = nansum(nansum(Cox_denitr(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(1000*hypsometry_Area(i-1));
        Hypsometry_oxid_rates_FeIII_molperm2(i-1) = nansum(nansum(Cox_FeIII(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(1000*hypsometry_Area(i-1));
        Hypsometry_oxid_rates_sulfate_molperm2(i-1) = nansum(nansum(Cox_sulfate(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(1000*hypsometry_Area(i-1));
        
        % SWI-fluxes for the hysometry in mumol cm-2 yr-1
        Hypsometry_SWI_flux_O2(i-1) = nansum(nansum(SWIflux_O2(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(hypsometry_Area(i-1));
        Hypsometry_SWI_flux_NO3(i-1) = nansum(nansum(SWIflux_NO3(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(hypsometry_Area(i-1));
        Hypsometry_SWI_flux_NH4(i-1) = nansum(nansum(SWIflux_NH4(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(hypsometry_Area(i-1));
        Hypsometry_SWI_flux_FeII(i-1) = nansum(nansum(SWIflux_FeII(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(hypsometry_Area(i-1));
        Hypsometry_SWI_flux_SO4(i-1) = nansum(nansum(SWIflux_SO4(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(hypsometry_Area(i-1));
        Hypsometry_SWI_flux_H2S(i-1) = nansum(nansum(SWIflux_H2S(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(hypsometry_Area(i-1));
        
        % Penetration depth for the hypsometry in cm
        Hypsometry_zox(i-1) = nansum(nansum(zox(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(hypsometry_Area(i-1));
        Hypsometry_zNO3(i-1) = nansum(nansum(zno3(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(hypsometry_Area(i-1));
        Hypsometry_zFeIII(i-1) = nansum(nansum(zfeIII(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(hypsometry_Area(i-1));
        Hypsometry_zSO4(i-1) = nansum(nansum(zso4(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(hypsometry_Area(i-1));
        
    end
    % Oxidation rates for the hypsometry in mumol cm-2 yr-1
    Hypsometry_oxid_rates_aerobic_mumolpercm2 = Hypsometry_oxid_rates_aerobic_molperm2*100^(-2)*10^6;
    Hypsometry_oxid_rates_denit_mumolpercm2 = Hypsometry_oxid_rates_denitr_molperm2*100^(-2)*10^6;
    Hypsometry_oxid_rates_FeIII_mumolpercm2 = Hypsometry_oxid_rates_FeIII_molperm2*100^(-2)*10^6;
    Hypsometry_oxid_rates_sulfate_mumolpercm2 = Hypsometry_oxid_rates_sulfate_molperm2*100^(-2)*10^6;
    
    
    % load Middelburg data
    filename = 'data/Middelburg96_GBCdataFig_1.xlsx';
    sheet = 'Rates_combined_DH';
    plot_data = false;
    
    Range_Aerobic_Atl = 'B6:C16';
    Aerobic_Atl = xlsread(filename,sheet,Range_Aerobic_Atl);
    Range_Aerobic_Pac = 'B25:C46';
    Aerobic_Pac = xlsread(filename,sheet,Range_Aerobic_Pac);
    Range_Aerobic_Ind = 'B51:C52';
    Aerobic_Ind = xlsread(filename,sheet,Range_Aerobic_Ind);
    
    Range_Denitr_Atl = 'F6:G21';
    Denitr_Atl = xlsread(filename,sheet,Range_Denitr_Atl);
    Range_Denitr_Pac = 'F25:G46';
    Denitr_Pac = xlsread(filename,sheet,Range_Denitr_Pac);
    Range_Denitr_Ind = 'F51:G52';
    Denitr_Ind = xlsread(filename,sheet,Range_Denitr_Ind);
    
    %% plot reaction rates
    set(0,'defaultLineLineWidth', 4)
    set(0,'DefaultAxesFontSize',26)
    
    fig_ox_rates = figure; %('Renderer', 'painters', 'Position', [10 10 2400 1200]);
    sgt = sgtitle('Oxidation rates (\mumol cm^{-2} yr^{-1})');
    set(gcf,'Position',[100 100 2000 1000])
    sgt.FontSize = 26;
    subplot(1,4,1)
    box on;
    hold on;
    plot(Hypsometry_oxid_rates_aerobic_mumolpercm2,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    if(plot_data)
        scatter(Aerobic_Atl(:,2),-Aerobic_Atl(:,1)/1000,90,'ks');
        scatter(Aerobic_Pac(:,2),-Aerobic_Pac(:,1)/1000,90,'k');
        scatter(Aerobic_Ind(:,2),-Aerobic_Ind(:,1)/1000,90,'k^');
    end
    xlim([0.0 20])
    ylim([-5.0 0.0])
    %    xlabel('Aerobic (\mumol cm^{-2} yr^{-1})')
    xlabel('Aerobic')
    ylabel('SFD  (km)');
    
    subplot(1,4,2)
    box on;
    hold on;
    plot(Hypsometry_oxid_rates_denit_mumolpercm2,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    if(plot_data)
    scatter(Denitr_Atl(:,2),-Denitr_Atl(:,1)/1000,90,'ks');
    scatter(Denitr_Pac(:,2),-Denitr_Pac(:,1)/1000,90,'k');
    scatter(Denitr_Ind(:,2),-Denitr_Ind(:,1)/1000,90,'k^');
    end
    xlim([0.0 20.0])
    ylim([-5.0 0.0])
    xlabel('Denitrification')
    ylabel('SFD  (km)');
    
    subplot(1,4,3)
    box on;
    hold on;
    plot(Hypsometry_oxid_rates_FeIII_mumolpercm2,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    xlim([0.0 5.0])
    ylim([-5.0 0.0])
    xlabel('FeIII reduction')
    ylabel('SFD  (km)');
    
    subplot(1,4,4)
    box on;
    hold on;
    plot(Hypsometry_oxid_rates_sulfate_mumolpercm2,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    xlim([0.0 200.0])
    ylim([-5.0 0.0])
    xlabel('Sulfate reduction')
    ylabel('SFD  (km)');
    
  	print(fig_ox_rates,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/Hypsometry_remin_rates_GAMMA50_' str_date '.eps']);    
    
    
    %% plot SWI-fluxes
    
    fig_SWI_fluxes = figure; %('Renderer', 'painters', 'Position', [10 10 2400 1200]);
    set(gcf,'Position',[100 100 1000 2000])
    sgt = sgtitle('SWI-fluxes (\mumol cm^{-2} yr^{-1})');
    sgt.FontSize = 26;
    subplot(2,3,1)
    box on;
    hold on;
    plot(Hypsometry_SWI_flux_O2,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    xlim([0.0 400])
    ylim([-5.0 0.0])
    %    xlabel('Aerobic (\mumol cm^{-2} yr^{-1})')
    xlabel('O_2')
    ylabel('SFD  (km)');
    
    subplot(2,3,2)
    box on;
    hold on;
    plot(Hypsometry_SWI_flux_NO3,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    xline(0,'k--','LineWidth',4);
    xlim([-3.5 3.5])
    ylim([-5.0 0.0])
    xlabel('NO_3')
    %    ylabel('SFD  (km)');
    
    subplot(2,3,3)
    box on;
    hold on;
    plot(Hypsometry_SWI_flux_SO4,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    xlim([0.0 10.0])
    ylim([-5.0 0.0])
    xlabel('SO_4')
    %    ylabel('SFD  (km)');
    
    subplot(2,3,4)
    box on;
    hold on;
    plot(Hypsometry_SWI_flux_NH4,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    xlim([-1.0 0.0])
    ylim([-5.0 0.0])
    xlabel('NH_4')
    ylabel('SFD  (km)');
    
    subplot(2,3,5)
    box on;
    hold on;
    plot(Hypsometry_SWI_flux_FeII,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    xlim([-0.1 0.0])
    ylim([-5.0 0.0])
    xlabel('Fe2+')
    %    ylabel('SFD  (km)');
    
    subplot(2,3,6)
    box on;
    hold on;
    plot(Hypsometry_SWI_flux_H2S,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    xlim([-10.0 0.0])
    ylim([-5.0 0.0])
    xlabel('H_2S')
    %    ylabel('SFD  (km)');
    
        print(fig_SWI_fluxes,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/Hypsometry_SWI_fluxes_GAMMA50_' str_date '.eps']);    
    
    
    %% plot penetration depths
    fig_penetration_depths = figure; %('Renderer', 'painters', 'Position', [10 10 1800 1200]);
    set(gcf,'Position',[100 100 450 900])
    box on;
    hold on;
    plot(Hypsometry_zox,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    plot(Hypsometry_zNO3,-hypsometry_meanDepth/1000,'rx--','MarkerFaceColor','r','MarkerSize',12);
%    plot(Hypsometry_zFeIII,-hypsometry_meanDepth/1000,'rx:','MarkerFaceColor','r','MarkerSize',12);
    %    plot(Hypsometry_zSO4,-hypsometry_meanDepth/1000,'k+-.','MarkerFaceColor','k','MarkerSize',12);
    %    xlim([0.0 100])
    ylim([-5.0 0.0])
    xlabel('Penetration depth (cm)','FontSize',24)
    ylabel('SFD  (km)','FontSize', 24);
    lgd = legend('zO_2','zNO_3');
%    lgd = legend('zO_2','zNO_3', 'zFeIII');
    lgd.FontSize = 20;
            print(fig_penetration_depths,'-depsc2', ['./plots/0504_ForKiel_Lee_a_PKfilled/Hypsometry_PenetrationDepths_GAMMA50_' str_date '.eps']);    

end


if(plot_BC)
    calc_res = 1;
    switch calc_res
        case 1
            load './data/O2_BW_WOA2018_hr_updated.mat'                  	% BW  O2 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
            load data/long.dat
            load data/lat.dat
            degree = '0.25degree';
            
            
        case 2
            load './data/BC_1degrees/O2_BW_WOA2018_01_updated.mat'                  	% BW  O2 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
            O2_BW_WOA2018_hr_updated = O2_BW_WOA2018_01_updated;
            degree = '1degree';
        case 3
            
    end
    
    fig_BWO2 = figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, O2_BW_WOA2018_hr_updated);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('BW O_2 (\mumol kg-1)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    %  caxis([-10.0 0.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_BWO2,'-depsc2', ['./plots/BW_O2_' degree '_' str_date '.eps']);
    
    
    
end

%%%%% Not needed

%% plot Fe2+ SWI-fluxes in Dale units
if(plot_Fe2_efflux)
    % in units of umol m-2 d-1 as in Dale ea. 2015, Fig. 3
    %     load output/Flux_Fe2_Dale_units_2m_1xTOC_gammaFe2_par_Fe3por0.85_SMIsinkDale_gamma0.95_101220.mat
    
    % units as Dale ea map Fig. 6a+c log(mmol m-2 yr-1)
    Flux_Fe2_logmmol = log10(Flux_Fe2_Dale_units*365*10^(-3));
    Flux_Fe2_logmmol(Flux_Fe2_logmmol==-Inf)=NaN;
    fig_Fe2_efflux = figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, Flux_Fe2_logmmol);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('Fe^{2+} SWI-flux (log(mmol m-2 yr-1))')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([-4 1.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_Fe2_efflux,'-depsc2', ['./plots/Fe2_efflux_' exp_name '_' str_date '_logmmolm2yr_until-4.eps']);
    
    % total efflux of Fe2+ in mol yr-1
    Global_Total_Fe2_efflux = nansum(nansum(Flux_Fe2_Dale_units.*dxdy))*10^(-6)*365;
    
    %
    %     fig_Fe2_efflux = figure;
    %     %set(gca,'FontSize',30)
    %     pcolor(long, lat, Flux_Fe2_Dale_units);
    %     hold on;
    %     % set(gca,'Color',lightGrey2)
    %     % set(gcf,'inverthardcopy','off');
    %     title('Fe^{2+} SWI-flux (\mumol m-2 d-1)')
    %     shading interp;
    %     %contour(long, lat, toc,'LineColor','k')
    %     colorbar ('horizontal')
    %     colormap(parula)
    %     caxis([0 1.0])
    %     xlabel('Longitude')
    %     ylabel('Latitude')
    %     print(fig_Fe2_efflux,'-depsc2', ['./plots/Fe2_efflux_' exp_name '_' str_date '_mumolm2day.eps']);
    %
    %     fig_Fe2_efflux = figure;
    %     %set(gca,'FontSize',30)
    %     pcolor(long, lat, Flux_Fe2_Dale_units*365);
    %     hold on;
    %     % set(gca,'Color',lightGrey2)
    %     % set(gcf,'inverthardcopy','off');
    %     title('Fe^{2+} SWI-flux (\mumol m-2 yr-1)')
    %     shading interp;
    %     %contour(long, lat, toc,'LineColor','k')
    %     colorbar ('horizontal')
    %     colormap(parula)
    %     caxis([0 200.0])
    %     xlabel('Longitude')
    %     ylabel('Latitude')
    %     print(fig_Fe2_efflux,'-depsc2', ['./plots/Fe2_efflux_' exp_name '_' str_date '_mumolm2yr.eps']);
    %
    
end

%% plot oxidation rates calculated from fluxes -- res.zTOC_RCM.calcCflx (saved in mmol m-2 day-1)
if(plot_oxidation_rates_fluxes)
    %    load output/Cox_rate_out_2m_1xTOC_gammaFe2_par_Fe3por0.85_SMIsinkDale_gamma0.95_101220.mat
    
    Cox_aerobic1 = -Cox_rate_out_flux.Cox_aerobic;
    Cox_denitr1=-Cox_rate_out_flux.Cox_denitr;
    Cox_FeIII1= -Cox_rate_out_flux.Cox_FeIII;
    Cox_sulfate1=-Cox_rate_out_flux.Cox_sulfate;
    Cox_total1=-Cox_rate_out_flux.Cox_Total;
    
    
    fig_OM_ox_total= figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, Cox_total1*1000);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('Total OM oxidation (\mumol m-2 d-1)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([0 200])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_OM_ox_total,'-depsc2', ['./plots/' exp_name '_OM_oxidation_Total_' str_date '_mummolm2day_until200.eps']);
    
    % Fraction of total
    fig_OM_ox_aerobic_frac= figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, Cox_aerobic1./Cox_total1*100);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('Fraction of aerobic degradation (in %)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([0 100])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_OM_ox_aerobic_frac,'-depsc2', ['./plots/' exp_name '_OM_oxidation_aerobic_frac_' str_date '_percent.eps']);
    
    
    fig_OM_ox_denit_frac= figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, Cox_denitr1./Cox_total1*100);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('Fraction of denitrification (in %)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([0 50])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_OM_ox_denit_frac,'-depsc2', ['./plots/' exp_name '_OM_oxidation_denit_frac_' str_date '_percent50.eps']);
    
    fig_OM_ox_Fe_red_frac= figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, Cox_FeIII1./Cox_total1*100);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('Fraction of Fe-reduction (in %)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([0 50])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_OM_ox_Fe_red_frac,'-depsc2', ['./plots/' exp_name '_OM_oxidation_Fe_red_frac_' str_date '_percent50.eps']);
    
    fig_OM_ox_sulfate_red_frac= figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, Cox_sulfate1./Cox_total1*100);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('Fraction of sulfate-reduction (in %)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([0 100])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_OM_ox_sulfate_red_frac,'-depsc2', ['./plots/' exp_name '_OM_oxidation_Sulfate_red_frac_' str_date '_percent.eps']);
    
end

% % end