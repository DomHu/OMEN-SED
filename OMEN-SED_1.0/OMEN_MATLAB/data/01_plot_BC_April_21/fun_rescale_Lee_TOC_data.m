% load original Lee et al. (2019) data in 5 × 5‐arc minute resolution ( = 1/12 x 1/12 degrees)
% and generate 1/4 x 1/4 ; 1 x 1 and 2 x 2 degree resolution; -- select the
% case for the resolution forther down in the code; i.e. calc_res = 
% NOTE: need to take different grid-cell size into account; 
% also check that at least 1/3% of grid-cells are ocean when calculating mean 
% otherwise define as land (i.e. set NaN)

load_and_plot_Lee_data = true;

path(path,'/home/domhu/Documents/MATLAB/M_Map');

if(load_and_plot_Lee_data)
    load('./data/DatasetS1_Lee_ea2019.txt');
    % reshape into matrix
    Lee_TOC=DatasetS1_Lee_ea2019(:,3);
    Lee_TOC_matrix = reshape(Lee_TOC,4320,2160);
    Lee_TOC_matrix = Lee_TOC_matrix';
    Lee_long= DatasetS1_Lee_ea2019(1:4320,1);
    Lee_lat= DatasetS1_Lee_ea2019(1:4320:end,2);

    % Robinson projection
    fig_toc_Lee = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    toc_levels = [0:0.1:5];
    [C,h] = m_contourf(Lee_long, Lee_lat, Lee_TOC_matrix,toc_levels);
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
 
end
    

%% Generate lower degree resolution
clear toc_sub dxdy_submatrix  Lee_toc_lr Lee_toc_lr_weighted

calc_res = 3;  % 1: 1/4; 2: 1 degree; 3: 2 degree resolution

% first calculate the area of the different grid-cells; needed to calculate the
% weighted average of TOC
dy = lat_to_m(1/12,Lee_lat);
dx = lon_to_m(1/12, Lee_lat);
dxdy   = dx.*dy;    % m^2 per grid-cell

dxdy_matrix = ones(size(Lee_TOC_matrix));
dxdy_matrix = dxdy_matrix.*dxdy;
% make sure dxxy_matrix=NaN where toc=NaN:
where_NaN = isnan(Lee_TOC_matrix);
dxdy_matrix(where_NaN)= NaN;

   
switch calc_res
    case 1
    	n_lower = 3;   % how much lower? 3 times: leads 1/4 degree resolution
    case 2
        n_lower = 12;   % how much lower? 12 times: leads 1degree resolution
    case 3
        n_lower = 24;   % how much lower? 24 times: leads 2 degree resolution
end
n_lower_help = n_lower - 1;    % for steps in matrix

[n,m]=size(Lee_TOC_matrix);
x_new = 0;
y_new = 0;

for x=1:n_lower:n-n_lower_help  % latitudes
    x_new=x_new+1;
    y_new = 0;
   for  y=1:n_lower:m-n_lower_help      % longitudes
      y_new = y_new+1;
      % toc submatrix
      toc_sub = Lee_TOC_matrix(x:x+n_lower_help, y:y+n_lower_help);               % n_lower x n_lower submatrix that will be combined
      % dxdy submatrix
      dxdy_submatrix = dxdy_matrix(x:x+n_lower_help, y:y+n_lower_help);               % n_lower x n_lower submatrix that will be combined
      % check if more than 66% is NaN:
      Num_NaNs = sum(sum(isnan(toc_sub)));
      NaN_max = 0.66*n_lower*n_lower;
      if(Num_NaNs > NaN_max)
    	Lee_toc_lr(x_new, y_new) = NaN;   % not enough cells with data -> NaN
        Lee_toc_lr_weighted(x_new, y_new) = NaN;   % not enough cells with data -> NaN

      else             
      	toc_nanmean =  nanmean(nanmean(toc_sub));     % calculate mean of toc
        toc_nanmean_weighted = nansum(nansum(toc_sub.*dxdy_submatrix))/nansum(nansum(dxdy_submatrix));
     	Lee_toc_lr(x_new, y_new) =toc_nanmean;            % save in new matrix
        Lee_toc_lr_weighted(x_new, y_new) =toc_nanmean_weighted;    % USE this as it is mean accounting for the grid-areas!
      end
      
   end
end

switch calc_res
    case 1
    save('./BC_Quarterdegree/Lee_toc_lr.mat' , 'Lee_toc_lr')      
    save('./BC_Quarterdegree/Lee_toc_lr_weighted.mat' , 'Lee_toc_lr_weighted')      
    lat_lr = (-89.875:0.25:89.875)';
    long_lr = (-179.875:0.25:179.875)';     
    save('./BC_Quarterdegree/lat_lr.mat' , 'lat_lr')      
    save('./BC_Quarterdegree/long_lr.mat' , 'long_lr')      
    

    case 2
    save('./BC_1degree/Lee_toc_lr.mat' , 'Lee_toc_lr')      
    save('./BC_1degree/Lee_toc_lr_weighted.mat' , 'Lee_toc_lr_weighted')      
    lat_lr = (-89.5:1.0:89.5)';
    long_lr = (-179.5:1.0:179.5)';
    save('./BC_1degree/lat_lr.mat' , 'lat_lr')      
    save('./BC_1degree/long_lr.mat' , 'long_lr')      

    case 3
    save('./BC_2degree/Lee_toc_lr.mat' , 'Lee_toc_lr')      
    save('./BC_2degree/Lee_toc_lr_weighted.mat' , 'Lee_toc_lr_weighted')      
    lat_lr = (-89.0:2.0:89.0)';
    long_lr = (-179.0:2.0:179.0)';
    save('./BC_2degree/lat_lr.mat' , 'lat_lr')      
    save('./BC_2degree/long_lr.mat' , 'long_lr')      

end

% plot lower resolution Lee TOC --  Robinson projection



    fig_toc_Lee_lr = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    toc_levels = [0:0.1:5];
    [C,h] = m_contourf(long_lr, lat_lr, Lee_toc_lr,toc_levels);
    set(h,'LineColor','none')
switch calc_res
    case 1
    title('Lee SWI-TOC 1/4 degree (wt%)')    
    case 2
    title('Lee SWI-TOC 1 degree (wt%)')
    case 3
        title('Lee SWI-TOC 2 degree (wt%)')
end
m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 3.0])
    xlabel('Longitude')
    ylabel('Latitude')
switch calc_res
    case 1
    print(fig_toc_Lee_lr,'-depsc2', ['TOC_input_Lee_Quartdegree.eps']);    
    case 2
    print(fig_toc_Lee_lr,'-depsc2', ['TOC_input_Lee_1degree.eps']);
    case 3
      print(fig_toc_Lee_lr,'-depsc2', ['TOC_input_Lee_2degree.eps']);  
end
    
    fig_toc_Lee_1degree_weighted = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    % set(gca,'FontSize',30)
    toc_levels = [0:0.1:5];
    [C,h] = m_contourf(long_lr, lat_lr, Lee_toc_lr_weighted,toc_levels);
    set(h,'LineColor','none')
switch calc_res
     case 1
    title('Lee SWI-TOC 1/4 degree weighted (wt%)')   
    case 2
    title('Lee SWI-TOC 1 degree weighted (wt%)')
    case 3
        title('Lee SWI-TOC 2 degree weighted (wt%)')
end
m_coast('linewidth',1,'color','k'); %('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    % hold on;
    % shading interp;
    colorbar ('horizontal')
    colormap(parula)
    caxis([0.0 3.0])
    xlabel('Longitude')
    ylabel('Latitude')
switch calc_res    
    case 1
    print(fig_toc_Lee_1degree_weighted,'-depsc2', ['TOC_input_Lee_Quartdegree_weighted.eps']);
    case 2
    print(fig_toc_Lee_1degree_weighted,'-depsc2', ['TOC_input_Lee_1degree_weighted.eps']);
    case 3
    print(fig_toc_Lee_1degree_weighted,'-depsc2', ['TOC_input_Lee_2degree_weighted.eps']);
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
