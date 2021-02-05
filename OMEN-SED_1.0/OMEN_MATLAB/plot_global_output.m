% function plot_global_output(high_res)
%% plots the results of the global analysis for OM degradation pathways
% high_res:     true    -   1/4 degree resolution
%               false   -   1 degree resolution

flags.calc_TOC_burial = false;
flags.calc_Regional_OM_oxidation = false;
flags.calc_Total_OM_oxidation = false;    % Global total in mol/yr -- taking size of grid-cells into account!
flags.plot_oxidation_rates = true;
flags.plot_penetration = true;
flags.plot_SWI_fluxes = true;
flags.calc_plot_hyp_remin_rates = true;

flags.plot_Fe2_efflux = false;
flags.plot_oxidation_rates_fluxes = false;  %% not needed - is the same as calculate rates with calcReac (flags.plot_oxidation_rates)

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
    load './data/BC_1degrees/water_depth_updated.mat'                           % Seafloor depth [m]
    water_depth_updated = -water_depth_updated;

    
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



if(flags.calc_TOC_burial)
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


if(flags.calc_Regional_OM_oxidation)
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

if(flags.calc_Total_OM_oxidation)
        
%    load output/dxdy_2m_1xTOC_gammaFe2_par_Fe3por0.85_SMIsinkDale_gamma0.95_101220.mat
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
if(flags.plot_oxidation_rates)
    % oxidation rates in mmol m-2 year!
%     load output/Cox_rate_out_2m_1xTOC_gammaFe2_par_Fe3por0.85_SMIsinkDale_gamma0.95_101220.mat
    
    Cox_aerobic = Cox_rate_out.Cox_aerobic;
    Cox_denitr=Cox_rate_out.Cox_denitr;
    Cox_FeIII= Cox_rate_out.Cox_FeIII;
    Cox_sulfate=Cox_rate_out.Cox_sulfate;
    Cox_total=Cox_rate_out.Cox_total;
    Total_rate_sum=Cox_rate_out.Total_rate_sum;
    
    
    fig_OM_ox_total= figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, Cox_total);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('Total OM oxidation (mmol m-2 yr-1)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([0 200])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_OM_ox_total,'-depsc2', ['./plots/' exp_name '_OM_oxidation_Total_' str_date '_mmmolm2yr_until200.eps']);
        
    % Fraction of total
	fig_OM_ox_aerobic_frac= figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, Cox_aerobic./Cox_total*100);
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
    pcolor(long, lat, Cox_denitr./Cox_total*100);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('Fraction of denitrification (in %)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([0 20])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_OM_ox_denit_frac,'-depsc2', ['./plots/' exp_name '_OM_oxidation_denit_frac_' str_date '_percent50.eps']);
        
	fig_OM_ox_Fe_red_frac= figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, Cox_FeIII./Cox_total*100);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('Fraction of Fe-reduction (in %)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([0 20])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_OM_ox_Fe_red_frac,'-depsc2', ['./plots/' exp_name '_OM_oxidation_Fe_red_frac_' str_date '_percent50.eps']);
        
	fig_OM_ox_sulfate_red_frac= figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, Cox_sulfate./Cox_total*100);
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




%% plot TEA penetration depth
if(flags.plot_penetration)
%    load output/Penetration_out_2m_1xTOC_gammaFe2_par_Fe3por0.85_SMIsinkDale_gamma0.95_101220.mat
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
    
    
    fig_SWI_BCond = figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, zox);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('O_2 penetration depth (cm)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    %caxis([-5000 0.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_SWI_BCond,'-depsc2', ['./plots/zox_' exp_name '_' str_date '.eps']);
    
    fig_SWI_BCond_no3 = figure;
    pcolor(long, lat, zno3);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('NO_3 penetration depth (cm)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([0 zinf])
    % pcolor(long, lat, zno3_inf);
    % shading interp;
    % colormap(gray)
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_SWI_BCond_no3,'-depsc2', ['./plots/zno3_' exp_name '_' str_date '.eps']);
    
    %
    
    fig_SWI_BCond_fe3 = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, zfeIII);
    title('Fe3^+ penetration depth (cm)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([0 zinf])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_SWI_BCond_fe3,'-depsc2', ['./plots/zfe3_' exp_name '_' str_date '.eps']);
    
    
    fig_SWI_BCond_so4 = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, zso4);
    title('SO_4 penetration depth (cm)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([0 zinf])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_SWI_BCond_so4,'-depsc2', ['./plots/zso4_' exp_name '_' str_date '.eps']);
    
end

%% plot calculated SWI-fluxes in mumol cm-2 yr-1
if(flags.plot_SWI_fluxes)
    
  	fig_SWI_fluxO2 = figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, SWI_fluxes.flxswiO2*10^6);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('SWI-flux O_2 (\mumol cm-2 yr-1)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([-50 0.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_SWI_fluxO2,'-depsc2', ['./plots/SWI_flux_O2_' exp_name '_' str_date '.eps']);

       
  	fig_SWI_fluxNO3 = figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, SWI_fluxes.flxswiNO3*10^6);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('SWI-flux NO_3 (\mumol cm-2 yr-1)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([-5 5.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_SWI_fluxNO3,'-depsc2', ['./plots/SWI_flux_NO3_' exp_name '_' str_date '.eps']);
        
  	fig_SWI_fluxNH4 = figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, SWI_fluxes.flxswiNH4*10^6);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('SWI-flux NH_4 (\mumol cm-2 yr-1)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 1.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_SWI_fluxNH4,'-depsc2', ['./plots/SWI_flux_NH4_' exp_name '_' str_date '.eps']);
 
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
    print(fig_Fe2_efflux,'-depsc2', ['./plots/SWI_flux_Fe2_' exp_name '_' str_date '_logmmolm2yr_until-4.eps']);

    
   	fig_SWI_fluxSO4 = figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, SWI_fluxes.flxswiSO4*10^6);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('SWI-flux SO_4 (\mumol cm-2 yr-1)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([-10.0 0.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_SWI_fluxSO4,'-depsc2', ['./plots/SWI_flux_SO4_' exp_name '_' str_date '.eps']);
     
   	fig_SWI_fluxH2S = figure;
    %set(gca,'FontSize',30)
    pcolor(long, lat, SWI_fluxes.flxswiH2S*10^6);
    hold on;
    % set(gca,'Color',lightGrey2)
    % set(gcf,'inverthardcopy','off');
    title('SWI-flux H_2S (\mumol cm-2 yr-1)')
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 10.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_SWI_fluxH2S,'-depsc2', ['./plots/SWI_flux_H2S_' exp_name '_' str_date '.eps']);
 
end


if(flags.calc_plot_hyp_remin_rates)
%% plot global hypsometry of remineralization rates (compare Fig. 5a-e, Thullner ea. 2009)
 
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
    
    % Oxidation rates for the hysometry in mol m-2 yr-1
    for i = 2:length(hypsometry_depths)
    Hypsometry_oxid_rates_aerobic_molperm2(i-1) = nansum(nansum(Cox_aerobic(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(1000*hypsometry_Area(i-1));
    Hypsometry_oxid_rates_denitr_molperm2(i-1) = nansum(nansum(Cox_denitr(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(1000*hypsometry_Area(i-1));
    Hypsometry_oxid_rates_FeIII_molperm2(i-1) = nansum(nansum(Cox_FeIII(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(1000*hypsometry_Area(i-1));
    Hypsometry_oxid_rates_sulfate_molperm2(i-1) = nansum(nansum(Cox_sulfate(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i)).*dxdy(water_depth_updated>=hypsometry_depths(i-1) & water_depth_updated<hypsometry_depths(i))))/(1000*hypsometry_Area(i-1));
    end
    % Oxidation rates for the hysometry in mumol cm-2 yr-1
    Hypsometry_oxid_rates_aerobic_mumolpercm2 = Hypsometry_oxid_rates_aerobic_molperm2*100^(-2)*10^6;
    Hypsometry_oxid_rates_denit_mumolpercm2 = Hypsometry_oxid_rates_denitr_molperm2*100^(-2)*10^6;
    Hypsometry_oxid_rates_FeIII_mumolpercm2 = Hypsometry_oxid_rates_FeIII_molperm2*100^(-2)*10^6;
    Hypsometry_oxid_rates_sulfate_mumolpercm2 = Hypsometry_oxid_rates_sulfate_molperm2*100^(-2)*10^6;

    
    % load Middelburg data
    filename = 'data/Middelburg96_GBCdataFig_1.xlsx';
    sheet = 'Rates_combined_DH';
    
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
        
    % plot reaction rates
    set(0,'defaultLineLineWidth', 4)
    set(0,'DefaultAxesFontSize',16)

    fig_ox_rates = figure('Renderer', 'painters', 'Position', [10 10 1800 1200]);
    subplot(1,4,1)
    box on;
    hold on;
    plot(Hypsometry_oxid_rates_aerobic_mumolpercm2,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    scatter(Aerobic_Atl(:,2),-Aerobic_Atl(:,1)/1000,90,'ks'); 
    scatter(Aerobic_Pac(:,2),-Aerobic_Pac(:,1)/1000,90,'k'); 
    scatter(Aerobic_Ind(:,2),-Aerobic_Ind(:,1)/1000,90,'k^'); 
    xlim([0.0 150])
    ylim([-5.0 0.0])
    xlabel('Aerobic (\mumol cm^{-2} yr^{-1})')
    ylabel('SFD  (km)');
    
    subplot(1,4,2)
    box on;
    hold on;
    plot(Hypsometry_oxid_rates_denit_mumolpercm2,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    scatter(Denitr_Atl(:,2),-Denitr_Atl(:,1)/1000,90,'ks'); 
    scatter(Denitr_Pac(:,2),-Denitr_Pac(:,1)/1000,90,'k'); 
    scatter(Denitr_Ind(:,2),-Denitr_Ind(:,1)/1000,90,'k^'); 
    xlim([0.0 20.0])
    ylim([-5.0 0.0])
    xlabel('Denitrification (\mumol cm^{-2} yr^{-1})')
    ylabel('SFD  (km)');
    
    subplot(1,4,3)
    box on;
    hold on;
    plot(Hypsometry_oxid_rates_FeIII_mumolpercm2,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    xlim([0.0 5.0])
    ylim([-5.0 0.0])
    xlabel('FeIII reduction (\mumol cm^{-2} yr^{-1})')
    ylabel('SFD  (km)');

    subplot(1,4,4)
    box on;
    hold on;
    plot(Hypsometry_oxid_rates_sulfate_mumolpercm2,-hypsometry_meanDepth/1000,'ko-','MarkerFaceColor','k','MarkerSize',12);
    % xlim([-5.0 2.0])
    ylim([-5.0 0.0])
    xlabel('Sulfate reduction (\mumol cm^{-2} yr^{-1})')
    ylabel('SFD  (km)');
    
    
    print(fig_ox_rates, '-depsc2', ['./plots/' exp_name '_hypsometry_remin_rates_' str_date '.eps']);

end


%%%%% Not needed

%% plot Fe2+ SWI-fluxes in Dale units
if(flags.plot_Fe2_efflux)
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
if(flags.plot_oxidation_rates_fluxes)
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