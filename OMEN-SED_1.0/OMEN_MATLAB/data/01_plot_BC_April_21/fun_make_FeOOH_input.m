% Calculate the total Fe flux as fct. of Holocene sedimentation rate
% compare Dale et al. (2015) Table 2:
% FluxFe_total = frac_fe * w * (1-por)*rho/Aw   % in mol /(cm2 yr)

%    load sed_holo.mat   % in cm/yr  -- 1/4 degree resolution
% load the updated sedimentation rate:
load('./BC_1degree/sed_lr_weighted.mat');

w = sed_lr_weighted;  % in cm/yr
frac_fe = 0.05;     %  Fe content in average sedimentary rock (~5%) [Garrels and Mackenzie, 1971] which is similar to Fe content in red clays [Glasby, 2006].
por = 0.85;          % porosity of compacted sediment
rho = 2.5;          % dry sediment density, in g cm-3)
Aw = 55.8;          %  standard atomic weight of iron (55.8 g mol-1)

FluxFe_total_Burw_BW_01_por85 = frac_fe * w * (1-por)*rho/Aw *10^6 ;                %  in mumol/(cm2 yr)
FluxFe_total_Burw_Daleunits_BW_01_por85 = FluxFe_total_Burw_BW_01_por85 *100^2 /365;    % in mumol/(m2 d))

por = 0.74;          % avg value as in Wallmann et al. (2012)
FluxFe_total_Burw_BW_01_por75 = frac_fe * w * (1-por)*rho/Aw *10^6 ;                %  in mumol/(cm2 yr)
FluxFe_total_Burw_Daleunits_BW_01_por74 = FluxFe_total_Burw_BW_01_por75 *100^2 /365;    % in mumol/(m2 d))

% save the BW concentrations
save('./BC_1degree/FluxFe_total_Burw_BW.mat' , 'FluxFe_total_Burw_BW_01_por85' , 'FluxFe_total_Burw_BW_01_por75')
save('./BC_1degree/FluxFe_total_Burw_Daleunits_BW.mat' , 'FluxFe_total_Burw_Daleunits_BW_01_por85', 'FluxFe_total_Burw_Daleunits_BW_01_por74')

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
[C,h] = m_contourf(long, lat, FluxFe_total_Burw_Daleunits_BW_01_por85,levels);
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

