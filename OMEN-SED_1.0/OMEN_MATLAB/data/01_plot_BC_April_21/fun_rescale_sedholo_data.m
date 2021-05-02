% load holocene sedimentation rate from Bradley ea. 2020 in 1/4 × 1/4 resolution
% the original data had a line and row too much:
% erase delete line 721 (just NaNs) and delete column 1441 and replace
% column 1 with 1440 (in old version column 1 are just NaNs)
% Theen generate 1 x 1 and 2 x 2 degree resolution; 
% -- select the case for the resolution forther down in the code; i.e. calc_res = 
% NOTE: need to take different grid-cell size into account; 
% also check that at least 1/3 of grid-cells are ocean when calculating mean 
% otherwise define as land (i.e. set NaN)

plot_original_sedrate = false;


path(path,'/home/domhu/Documents/MATLAB/M_Map');


% load original sedimentation rate
% old lat/long is wrong; should have the average lat/long of every grid-cell
    load('./data/Bradley/sed_holo.mat');        

if(plot_original_sedrate)
    load('./data/Bradley/long.dat');
    load('./data/Bradley/lat.dat');
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
end

% clean original sed-rate data
% delete line 721
    sed_holo_clean = sed_holo(1:end-1,1:end-1);
    sed_holo_clean(:,1) = sed_holo_clean(:,end);
    
    % plot updated sedimentation rate
   	lat_quarter = (-89.875:0.25:89.875)';
    long_quarter = (-179.875:0.25:179.875)';     

   	% using Robinson projection
    fig_sedrate_updated = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    sed_levels = [0:0.005:0.13];
    [C,h] = m_contourf(long_quarter, lat_quarter, sed_holo_clean, sed_levels);    % 1/4 degree resolution
    set(h,'LineColor','none')
    title('Sedimentation rate clean 1/4^\circ (cm/yr)')
    m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 0.12])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_sedrate_updated,'-depsc2', ['Sed_Rate_clean_Quartdegree.eps']);
    save('./BC_Quarterdegree/sed_holo_clean.mat' , 'sed_holo_clean')      

    
%% Generate lower degree resolution
clear sed_lr sed_lr_weighted

calc_res = 2;  % 1: 1 degree; 2: 2 degree resolution

% first calculate the area of the different grid-cells; needed to calculate the
% weighted average of TOC
dy = lat_to_m(1/12,lat_quarter);
dx = lon_to_m(1/12, lat_quarter);
dxdy   = dx.*dy;    % m^2 per grid-cell

dxdy_matrix = ones(size(sed_holo_clean));
dxdy_matrix = dxdy_matrix.*dxdy;
% make sure dxxy_matrix=NaN where sed_holo_clean=NaN:
where_NaN = isnan(sed_holo_clean);
dxdy_matrix(where_NaN)= NaN;

switch calc_res
    case 1
    	n_lower = 4;   % how much lower? 4 times: leads 1 degree resolution
    case 2
        n_lower = 8;   % how much lower? 8 times: leads 2 degree resolution
end
n_lower_help = n_lower - 1;    % for steps in matrix

[n,m]=size(sed_holo_clean);
x_new = 0;
y_new = 0;

for x=1:n_lower:n-n_lower_help  % latitudes
    x_new=x_new+1;
    y_new = 0;
   for  y=1:n_lower:m-n_lower_help      % longitudes
      y_new = y_new+1;
      % toc submatrix
      sed_sub = sed_holo_clean(x:x+n_lower_help, y:y+n_lower_help);               % n_lower x n_lower submatrix that will be combined
      % dxdy submatrix
      dxdy_submatrix = dxdy_matrix(x:x+n_lower_help, y:y+n_lower_help);               % n_lower x n_lower submatrix that will be combined
      % check if more than 66% is NaN:
      Num_NaNs = sum(sum(isnan(sed_sub)));
      NaN_max = 0.66*n_lower*n_lower;
      if(Num_NaNs > NaN_max)
    	sed_lr(x_new, y_new) = NaN;   % not enough cells with data -> NaN
        sed_lr_weighted(x_new, y_new) = NaN;   % not enough cells with data -> NaN

      else             
      	sed_nanmean =  nanmean(nanmean(sed_sub));     % calculate mean of toc
        sed_nanmean_weighted = nansum(nansum(sed_sub.*dxdy_submatrix))/nansum(nansum(dxdy_submatrix));
     	sed_lr(x_new, y_new) =sed_nanmean;            % save in new matrix
        sed_lr_weighted(x_new, y_new) =sed_nanmean_weighted;
      end
      
   end
end

switch calc_res
    case 1
    save('./BC_1degree/sed_lr.mat' , 'sed_lr')      
    save('./BC_1degree/sed_lr_weighted.mat' , 'sed_lr_weighted')      
    lat_lr = (-89.5:1.0:89.5)';
    long_lr = (-179.5:1.0:179.5)';

    case 2
    save('./BC_2degree/sed_lr.mat' , 'sed_lr')      
    save('./BC_2degree/sed_lr_weighted.mat' , 'sed_lr_weighted')      
    lat_lr = (-89.0:2.0:89.0)';
    long_lr = (-179.0:2.0:179.0)';

end

% plot lower resolution Lee TOC --  Robinson projection

    fig_sed_lr = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    sed_levels = [0:0.005:0.13];
    [C,h] = m_contourf(long_lr, lat_lr, sed_lr, sed_levels);
    set(h,'LineColor','none')
switch calc_res
    case 1
    title('Sedimentation rate clean 1^\circ (cm/yr)')
    case 2
    title('Sedimentation rate clean 2^\circ (cm/yr)')
end
m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 0.12])
    xlabel('Longitude')
    ylabel('Latitude')
switch calc_res
    case 1
    print(fig_sed_lr,'-depsc2', ['Sed_Rate_clean_1degree.eps']);
    case 2
      print(fig_sed_lr,'-depsc2', ['Sed_Rate_clean_2degree.eps']);  
end
    
    fig_sed_lr_weighted = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    sed_levels = [0:0.005:0.13];
    [C,h] = m_contourf(long_lr, lat_lr, sed_lr_weighted,sed_levels);
    set(h,'LineColor','none')
switch calc_res
    case 1
    title('Sedimentation rate clean 1^\circ weighted (cm/yr)')
    case 2
    title('Sedimentation rate clean 2^\circ weighted (cm/yr)')
end
m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 0.12])
    xlabel('Longitude')
    ylabel('Latitude')
switch calc_res    
    case 1
    print(fig_sed_lr_weighted,'-depsc2', ['Sed_Rate_clean_1degree_weighted.eps']);
    case 2
    print(fig_sed_lr_weighted,'-depsc2', ['Sed_Rate_clean_2degree_weighted.eps']);
end

    
    
 

function dy = lat_to_m(dlat,alat)
%  dy = lat_to_m(dlat,alat)
% dy   = latitude difference in meters
% dlat = latitude difference in degrees
% alat = average latitude between the two fixes
% Reference: American Practical Navigator, Vol II, 1975 Edition, p 5 

rlat = alat * pi/180;
m = 111132.09 * ones(size(rlat)) - ...
    566.05 * cos(2 * rlat) + 1.2 * cos(4 * rlat);
dy = dlat .* m ;
end


function dx = lon_to_m(dlon, alat)
% dx = lon_to_m(dlon, alat)
% dx   = longitude difference in meters
% dlon = longitude difference in degrees
% alat = average latitude between the two fixes

rlat = alat * pi/180;
p = 111415.13 * cos(rlat) - 94.55 * cos(3 * rlat);
dx = dlon .* p;
end
