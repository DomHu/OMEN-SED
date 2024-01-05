% original file was used in ~/Documents/Projects/18_RECCAP2/OMEN/DOU_maps

path(path,'/home/domhu/Documents/MATLAB/M_Map');

plot_BC = true;
fill_up_missing_DOU = true;
plot_DOU_mask_coast_filled = true;

lat = (89.875:-0.25:-89.875);
long = (-179.875:0.25:179.875);


load('./DOU_mask.mat');

toc = load('../Markus_data_AprilMay23/surfOC_matrix_2023_new_v01.csv');

Area_mask = toc;
Area_mask(~isnan(toc))=1;

water_depth =	-load('../Markus_data_April22/bathymetry_matrix_new.csv');
water_depth = water_depth.*Area_mask;



if plot_BC
    fig_DOU = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:0.5:15];
    limits = [0 10];
    [C,h] = m_contourf(long, lat,  DOU_mask, levels);
    set(h,'LineColor','none')
    title('DOU (mmol m-2 d-1)') % -- Lee et al. (2019)')
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    xlabel('Longitude')
    ylabel('Latitude')
    
    
    fig_toc = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:0.1:5.0];
    limits = [0 3.0];
    [C,h] = m_contourf(long, lat,  toc, levels);
    set(h,'LineColor','none')
    title('TOC (wt%)') % -- Lee et al. (2019)')
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    xlabel('Longitude')
    ylabel('Latitude')
end

if fill_up_missing_DOU
    
    DOU_mask_coast_filled = DOU_mask;
    DOU_mask_coast_filled = DOU_mask_coast_filled.*Area_mask;

    % create defeault DOU related to grid depth
    DOU_predicted_depthbin = [20.0 8.0 3.1 1.4 0.8 0.6]; % in mmol O2 m-2 d-1 from Jorgensen et al. (2022) Table 4
  	DOU_griddepth = NaN(size(DOU_mask));

    depth_bins = [50 200 1000 2000 4000];
    DOU_griddepth(water_depth<=depth_bins(1)) = DOU_predicted_depthbin(1);
    DOU_griddepth(water_depth<=depth_bins(2) & water_depth>depth_bins(1)) = DOU_predicted_depthbin(2);
    DOU_griddepth(water_depth<=depth_bins(3) & water_depth>depth_bins(2)) = DOU_predicted_depthbin(3);
    DOU_griddepth(water_depth<=depth_bins(4) & water_depth>depth_bins(3)) = DOU_predicted_depthbin(4);
    DOU_griddepth(water_depth<=depth_bins(5) & water_depth>depth_bins(4)) = DOU_predicted_depthbin(5);
    DOU_griddepth(water_depth>depth_bins(5)) = DOU_predicted_depthbin(6);

    
    % matrix: true = DOU-value missing
    DOU_missing = (~isnan(toc) & isnan(DOU_mask));
    % 1st fill-up with mean of neighboring cells
	[m,n]=size(toc);
            
    % use nanmean of 4 neighboring cells, i.e. ~ 0.5° in each direction
    % if still NaN use water depth average
    for x = 3:m-2
        for y = 3:n-2
            if DOU_missing(x,y)==1
                DOU_mask_coast_filled(x,y)=nanmean(reshape(DOU_mask(x-2:x+2,y-2:y+2),1,[]));  % use here only the original DOU_mask values -- otherwise 
                if isnan(DOU_mask_coast_filled(x,y))
                    DOU_mask_coast_filled(x,y) = DOU_griddepth(x,y);
                end
            end
        end
    end
    
	DOU_missing = (~isnan(toc) & isnan(DOU_mask_coast_filled));
    DOU_mask_coast_filled(DOU_missing) = DOU_griddepth(DOU_missing);
    
	save('./DOU_mask_coast_filled.mat' , 'DOU_mask_coast_filled')

    if plot_DOU_mask_coast_filled
    
    fig_DOU_filled = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:0.5:15];
    limits = [0 10];
    [C,h] = m_contourf(long, lat,  DOU_mask_coast_filled, levels);
    set(h,'LineColor','none')
    title('DOU filled (mmol m-2 d-1)') % -- Lee et al. (2019)')
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    xlabel('Longitude')
    ylabel('Latitude')

  	print(fig_DOU_filled,'-dpdf', ['./01_Mask_DOU_mask_coast_filled.pdf']);    
   % exportgraphics(gcf,['./01_Mask_DOU_mask_coast_filled_export.pdf']);    % Since R2020TLAB


    end

end

















