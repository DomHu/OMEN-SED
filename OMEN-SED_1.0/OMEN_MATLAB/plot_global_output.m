% function plot_global_output(high_res)
%% plots the results of the global analysis for OM degradation pathways
% high_res:     true    -   1/4 degree resolution
%               false   -   1 degree resolution

calc_TOC_burial = true;
calc_Regional_OM_oxidation = false;
calc_Total_OM_oxidation = false;    % Global total in mol/yr -- taking size of grid-cells into account!
plot_oxidation_rates = false;
plot_penetration = false;
plot_SWI_fluxes = false;

plot_Fe2_efflux = false;
plot_oxidation_rates_fluxes = false;  %% not needed - is the same as calculate rates with calcReac (plot_oxidation_rates)

zinf = 50;

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
    Total_OM_burial_Pg_Global = 12*(Total_OM_burial_mol_Global)*10^(-15);               

    %% Depth regime burial values in Pg C yr-1
    
    Total_OM_burial_Pg_NearShore = nansum(nansum(TOC_burial_flux_m2(water_depth_updated<near_shore).*dxdy(water_depth_updated<near_shore)))*12*10^(-15);
    Total_OM_burial_Pg_Shelf = nansum(nansum(TOC_burial_flux_m2(water_depth_updated>=near_shore & water_depth_updated<shelf).*dxdy(water_depth_updated>=near_shore & water_depth_updated<shelf)))*12*10^(-15);
    Total_OM_burial_Pg_Slope = nansum(nansum(TOC_burial_flux_m2(water_depth_updated>=shelf & water_depth_updated<slope).*dxdy(water_depth_updated>=shelf & water_depth_updated<slope)))*12*10^(-15);
    Total_OM_burial_Pg_Plain = nansum(nansum(TOC_burial_flux_m2(water_depth_updated>=slope).*dxdy(water_depth_updated>=slope)))*12*10^(-15);

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
    Total_OM_oxid_Pg_Global = 12*(Total_OM_oxid_mmol_Global/1000)*10^(-15);               

    %% Total depth regime oxidation rates in mmol yr-1
    Total_OM_oxid_mmol_NearShore = nansum(nansum(Cox_total(water_depth_updated<near_shore).*dxdy(water_depth_updated<near_shore)));
   	Total_OM_oxid_mmol_Shelf = nansum(nansum(Cox_total(water_depth_updated>=near_shore & water_depth_updated<shelf).*dxdy(water_depth_updated>=near_shore & water_depth_updated<shelf)));
   	Total_OM_oxid_mmol_Slope = nansum(nansum(Cox_total(water_depth_updated>=shelf & water_depth_updated<slope).*dxdy(water_depth_updated>=shelf & water_depth_updated<slope)));
   	Total_OM_oxid_mmol_Plain = nansum(nansum(Cox_total(water_depth_updated>=slope).*dxdy(water_depth_updated>=slope)));

    % Total depth regime oxidation rates in Pg yr-1 (from mmol yr-1) 
    Total_OM_oxid_Pg_NearShore = 12*(Total_OM_oxid_mmol_NearShore/1000)*10^(-15);               
    Total_OM_oxid_Pg_Shelf = 12*(Total_OM_oxid_mmol_Shelf/1000)*10^(-15);               
    Total_OM_oxid_Pg_Slope = 12*(Total_OM_oxid_mmol_Slope/1000)*10^(-15);               
    Total_OM_oxid_Pg_Plain = 12*(Total_OM_oxid_mmol_Plain/1000)*10^(-15);              
    
    
    %% TEA specific depth regime oxidation rates in Pg C yr-1
    
    % Near shore
 	NearShore_aerobic_OM_oxid_Pg = nansum(nansum(Cox_aerobic(water_depth_updated<near_shore).*dxdy(water_depth_updated<near_shore)))*12/1000*10^(-15);
 	NearShore_denitr_OM_oxid_Pg = nansum(nansum(Cox_denitr(water_depth_updated<near_shore).*dxdy(water_depth_updated<near_shore)))*12/1000*10^(-15);
 	NearShore_FeIII_OM_oxid_Pg = nansum(nansum(Cox_FeIII(water_depth_updated<near_shore).*dxdy(water_depth_updated<near_shore)))*12/1000*10^(-15);
 	NearShore_sulfate_OM_oxid_Pg = nansum(nansum(Cox_sulfate(water_depth_updated<near_shore).*dxdy(water_depth_updated<near_shore)))*12/1000*10^(-15);

    Frac_NearShore_aerobic = NearShore_aerobic_OM_oxid_Pg/Total_OM_oxid_Pg_NearShore*100;
 	Frac_NearShore_denitr = NearShore_denitr_OM_oxid_Pg/Total_OM_oxid_Pg_NearShore*100;
 	Frac_NearShore_FeIII = NearShore_FeIII_OM_oxid_Pg/Total_OM_oxid_Pg_NearShore*100;
 	Frac_NearShore_sulfate = NearShore_sulfate_OM_oxid_Pg/Total_OM_oxid_Pg_NearShore*100;

    if(abs(NearShore_aerobic_OM_oxid_Pg+NearShore_denitr_OM_oxid_Pg+NearShore_FeIII_OM_oxid_Pg+NearShore_sulfate_OM_oxid_Pg - Total_OM_oxid_Pg_NearShore)>1e-6)
        error('NearShore ~= 1.0');
    end

    % Shelf
 	Shelf_aerobic_OM_oxid_Pg = nansum(nansum(Cox_aerobic(water_depth_updated>=near_shore & water_depth_updated<shelf).*dxdy(water_depth_updated>=near_shore & water_depth_updated<shelf)))*12/1000*10^(-15);
 	Shelf_denitr_OM_oxid_Pg = nansum(nansum(Cox_denitr(water_depth_updated>=near_shore & water_depth_updated<shelf).*dxdy(water_depth_updated>=near_shore & water_depth_updated<shelf)))*12/1000*10^(-15);
 	Shelf_FeIII_OM_oxid_Pg = nansum(nansum(Cox_FeIII(water_depth_updated>=near_shore & water_depth_updated<shelf).*dxdy(water_depth_updated>=near_shore & water_depth_updated<shelf)))*12/1000*10^(-15);
 	Shelf_sulfate_OM_oxid_Pg = nansum(nansum(Cox_sulfate(water_depth_updated>=near_shore & water_depth_updated<shelf).*dxdy(water_depth_updated>=near_shore & water_depth_updated<shelf)))*12/1000*10^(-15);

    Frac_Shelf_aerobic = Shelf_aerobic_OM_oxid_Pg/Total_OM_oxid_Pg_Shelf*100;
 	Frac_Shelf_denitr = Shelf_denitr_OM_oxid_Pg/Total_OM_oxid_Pg_Shelf*100;
 	Frac_Shelf_FeIII = Shelf_FeIII_OM_oxid_Pg/Total_OM_oxid_Pg_Shelf*100;
 	Frac_Shelf_sulfate = Shelf_sulfate_OM_oxid_Pg/Total_OM_oxid_Pg_Shelf*100;
        
    if(abs(Shelf_aerobic_OM_oxid_Pg+Shelf_denitr_OM_oxid_Pg+Shelf_FeIII_OM_oxid_Pg+Shelf_sulfate_OM_oxid_Pg - Total_OM_oxid_Pg_Shelf)>1e-6)
        error('Shelf ~= 1.0');
    end

    % Slope
 	Slope_aerobic_OM_oxid_Pg = nansum(nansum(Cox_aerobic(water_depth_updated>=shelf & water_depth_updated<slope).*dxdy(water_depth_updated>=shelf & water_depth_updated<slope)))*12/1000*10^(-15);
 	Slope_denitr_OM_oxid_Pg = nansum(nansum(Cox_denitr(water_depth_updated>=shelf & water_depth_updated<slope).*dxdy(water_depth_updated>=shelf & water_depth_updated<slope)))*12/1000*10^(-15);
 	Slope_FeIII_OM_oxid_Pg = nansum(nansum(Cox_FeIII(water_depth_updated>=shelf & water_depth_updated<slope).*dxdy(water_depth_updated>=shelf & water_depth_updated<slope)))*12/1000*10^(-15);
 	Slope_sulfate_OM_oxid_Pg = nansum(nansum(Cox_sulfate(water_depth_updated>=shelf & water_depth_updated<slope).*dxdy(water_depth_updated>=shelf & water_depth_updated<slope)))*12/1000*10^(-15);

    Frac_Slope_aerobic = Slope_aerobic_OM_oxid_Pg/Total_OM_oxid_Pg_Slope*100;
 	Frac_Slope_denitr = Slope_denitr_OM_oxid_Pg/Total_OM_oxid_Pg_Slope*100;
 	Frac_Slope_FeIII = Slope_FeIII_OM_oxid_Pg/Total_OM_oxid_Pg_Slope*100;
 	Frac_Slope_sulfate = Slope_sulfate_OM_oxid_Pg/Total_OM_oxid_Pg_Slope*100;
    
    if(abs(Slope_aerobic_OM_oxid_Pg+Slope_denitr_OM_oxid_Pg+Slope_FeIII_OM_oxid_Pg+Slope_sulfate_OM_oxid_Pg - Total_OM_oxid_Pg_Slope)>1e-6)
        error('Slope ~= 1.0');
    end

    % Plain
 	Plain_aerobic_OM_oxid_Pg = nansum(nansum(Cox_aerobic(water_depth_updated>=slope).*dxdy(water_depth_updated>=slope)))*12/1000*10^(-15);
 	Plain_denitr_OM_oxid_Pg = nansum(nansum(Cox_denitr(water_depth_updated>=slope).*dxdy(water_depth_updated>=slope)))*12/1000*10^(-15);
 	Plain_FeIII_OM_oxid_Pg = nansum(nansum(Cox_FeIII(water_depth_updated>=slope).*dxdy(water_depth_updated>=slope)))*12/1000*10^(-15);
 	Plain_sulfate_OM_oxid_Pg = nansum(nansum(Cox_sulfate(water_depth_updated>=slope).*dxdy(water_depth_updated>=slope)))*12/1000*10^(-15);

    Frac_Plain_aerobic = Plain_aerobic_OM_oxid_Pg/Total_OM_oxid_Pg_Plain*100;
 	Frac_Plain_denitr = Plain_denitr_OM_oxid_Pg/Total_OM_oxid_Pg_Plain*100;
 	Frac_Plain_FeIII = Plain_FeIII_OM_oxid_Pg/Total_OM_oxid_Pg_Plain*100;
 	Frac_Plain_sulfate = Plain_sulfate_OM_oxid_Pg/Total_OM_oxid_Pg_Plain*100;
    
    if(abs(Plain_aerobic_OM_oxid_Pg+Plain_denitr_OM_oxid_Pg+Plain_FeIII_OM_oxid_Pg+Plain_sulfate_OM_oxid_Pg - Total_OM_oxid_Pg_Plain)>1e-6)
        error('Plain ~= 1.0');
    end
    
    %% calculate area of depth regimes
    Area_NearShore = nansum(nansum(dxdy(water_depth_updated<near_shore)));
    Area_Shelf = nansum(nansum(dxdy(water_depth_updated>=near_shore & water_depth_updated<shelf)));
    Area_Slope = nansum(nansum(dxdy(water_depth_updated>=shelf & water_depth_updated<slope)));
    Area_Plain = nansum(nansum(dxdy(water_depth_updated>=slope)));

end

if(calc_Total_OM_oxidation)
        
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
    Global_Total_OM_oxid = nansum(nansum(Cox_total.*dxdy));
    Global_Total_OM_oxid_Tmolyr = (Global_Total_OM_oxid/1000)/10^(12);
	Global_aerobic_OM_oxid = nansum(nansum(Cox_aerobic.*dxdy));
	Global_denitr_OM_oxid = nansum(nansum(Cox_denitr.*dxdy));
	Global_FeIII_OM_oxid = nansum(nansum(Cox_FeIII.*dxdy));
	Global_sulfate_OM_oxid = nansum(nansum(Cox_sulfate.*dxdy));

    Global_aerobic_Frac = Global_aerobic_OM_oxid/Global_Total_OM_oxid*100;
   	Global_denit_Frac = Global_denitr_OM_oxid/Global_Total_OM_oxid*100;
    Global_Fe3_Frac = Global_FeIII_OM_oxid/Global_Total_OM_oxid*100;
    Global_sulfate_Frac = Global_sulfate_OM_oxid/Global_Total_OM_oxid*100;
    
    % oxidation rate in g yr-1 (from mmol yr-1) 
    Global_Total_OM_oxid_gram = 12*Global_Total_OM_oxid/1000;               
     
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
    caxis([0 50])
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
    caxis([0 50])
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
if(plot_penetration)
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
if(plot_SWI_fluxes)
    
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