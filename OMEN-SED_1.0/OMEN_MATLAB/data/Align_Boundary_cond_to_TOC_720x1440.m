%% Make sure every grid-cell with TOC-values also has a value for the BC
%% Here using TOC for the RECCAP2 manuscript - hence, only TOC for the margins
% using the higest resolution in 1/4 degrees
% use mean of neighboring cells - might need to run it multiple times to
% assure all grid-cells are filled with a value
path(path,'/home/domhu/Documents/MATLAB/M_Map');


plot_BC_in = false;

align_Bcs = true;

str_date = [datestr(date,7), datestr(date,5), datestr(date,11)];

lat = (89.875:-0.25:-89.875);
long = (-179.875:0.25:179.875);

toc = load('./BC_calc_a_from_Jorgensen/surfOC_matrix_2023_new_v01.csv');

toc_values = ~isnan(toc);

BC_Tmp_in = flipud(struct2array(load('./WOA2018/Tmp_BW_WOA2018_hr.mat')));
BC_O2_in = flipud(struct2array(load('./WOA2018/O2_BW_WOA2018_hr.mat')));
BC_NO3_in = flipud(struct2array(load('./WOA2018/NO3_BW_WOA2018_hr.mat')));

% Boundary conditions just for our margin areas
BC_Tmp_Margin = BC_Tmp_in;
BC_Tmp_Margin(isnan(toc)) = NaN;

BC_O2_Margin = BC_O2_in;
BC_O2_Margin(isnan(toc)) = NaN;

BC_NO3_Margin = BC_NO3_in;
BC_NO3_Margin(isnan(toc)) = NaN;

fig(1) = figure;
m_proj('Robinson','longitudes',[-180 179.99], ...
    'latitudes',[-90 90]);
hold on;
levels = [0:0.1:5.0];
limits = [0 3.0];
m_coast('linewidth',1,'color','k');
m_coast('patch',[.7 .7 .7]);
[C,h] = m_contourf(long, lat,  toc, levels);
set(h,'LineColor','none')
title('TOC (wt%) ') % -- Lee et al. (2019)')
m_grid('fontsize',8);
colorbar ('horizontal')
colormap(parula)
caxis(limits)
xlabel('Longitude')
ylabel('Latitude')


if plot_BC_in
    
%     fig(2) = figure;
%     m_proj('Robinson','longitudes',[-180 179.99], ...
%         'latitudes',[-90 90]);
%     hold on;
%     levels = [0:0.5:21.0];
%     limits = [0 20.0];
%     m_coast('linewidth',1,'color','k');
%     m_coast('patch',[.7 .7 .7]);
%     [C,h] = m_contourf(long, lat,  BC_Tmp_in, levels);
%     set(h,'LineColor','none')
%     title('Seafloor temperature (Celsius) ') % -- Lee et al. (2019)')
%     m_grid('fontsize',8);
%     colorbar ('horizontal')
%     colormap(parula)
%     caxis(limits)
%     xlabel('Longitude')
%     ylabel('Latitude')
    
    fig(21) = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:0.5:21.0];
    limits = [0 20.0];
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    [C,h] = m_contourf(long, lat,  BC_Tmp_Margin, levels);
    set(h,'LineColor','none')
    title('Seafloor temperature (Celsius) in') % -- Lee et al. (2019)')
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    xlabel('Longitude')
    ylabel('Latitude')
    
    
    fig(3) = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:10:310.0];
    limits = [0 300.0];
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    [C,h] = m_contourf(long, lat,  BC_O2_Margin, levels);
    set(h,'LineColor','none')
    title('Seafloor Oxygen (\muM) in ') % -- Lee et al. (2019)')
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    xlabel('Longitude')
    ylabel('Latitude')
    
    
    
    fig(4) = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:2:52.0];
    limits = [0 50.0];
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    [C,h] = m_contourf(long, lat,  BC_NO3_Margin, levels);
    set(h,'LineColor','none')
    title('Seafloor Nitrate (\muM) in') % -- Lee et al. (2019)')
    m_grid('fontsize',8);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    xlabel('Longitude')
    ylabel('Latitude')
    
end



toc_notNaN = (~isnan(toc));
size_toc = sum(sum(toc_notNaN));
[m,n]=size(toc);


% loop through the different boundary conditions (i.e. 8 different ones)
if(align_Bcs)
    
    for bc=1:3
        switch bc
            case 1  % Temperature
                BCond_BW_hr = BC_Tmp_Margin;
                BC_global = BC_Tmp_in; % to calculate local averages
               	txt = 'Temperature';
            case 2  % O2
                BCond_BW_hr = BC_O2_Margin;
                BC_global = BC_O2_in; % to calculate local averages
               	txt = 'Oxygen';
            case 3  % NO3
                BCond_BW_hr = BC_NO3_Margin;
                BC_global = BC_NO3_in; % to calculate local averages
               	txt = 'Nitrate';
        end
        
        
        % use nanmean of 4 neighboring cells, i.e. ~ 0.25° in each direction
        % if still NaN use water depth average
      	BC_missing = (~isnan(toc) & isnan(BCond_BW_hr));
        Sum_BC_mising_new = sum(sum(BC_missing)); 
        Sum_BC_mising_old = 0;
        iter = 0;
        
        while Sum_BC_mising_new ~= Sum_BC_mising_old
            for x = 3:m-2
                for y = 3:n-2
                    if BC_missing(x,y)==1
                        BCond_BW_hr(x,y)=nanmean(reshape(BC_global(x-2:x+2,y-2:y+2),1,[]));  % use here only the original DOU_mask values -- otherwise
                        %                     if isnan(DOU_mask_coast_filled(x,y))
                        %                         DOU_mask_coast_filled(x,y) = DOU_griddepth(x,y);
                        %                     end
                    end
                end
            end
            Sum_BC_mising_old = Sum_BC_mising_new;
            BC_missing = (~isnan(toc) & isnan(BCond_BW_hr));
            Sum_BC_mising_new = sum(sum(BC_missing));            
        end
        % Cells with still NaN use global mean
        Mean_value = nanmean(nanmean(BCond_BW_hr));
        BCond_BW_hr(BC_missing)= Mean_value;
        missing_perc = 100*Sum_BC_mising_new/size_toc;
      	fprintf([txt, ' %i missing (or %4.2f percent), get mean value of %4.2f\n'],  Sum_BC_mising_new, missing_perc, Mean_value);

        
        switch bc
            case 1
                BC_Tmp_Margin_updated = BCond_BW_hr;
                save('BC_Tmp_Margin_updated.mat' , 'BC_Tmp_Margin_updated')
                save('./BC_calc_a_from_Jorgensen/BC_Tmp_Margin_updated.mat' , 'BC_Tmp_Margin_updated')
                txt_file = 'Temp_WOA_';
                txt_title = 'Temperature at the SWI (\circC)';
                 levels = [0:0.5:21.0];
                 limits = [0 20.0];
                
            case 2  % O2
                BC_O2_Margin_updated = BCond_BW_hr;
                save('BC_O2_Margin_updated.mat' , 'BC_O2_Margin_updated')
                save('./BC_calc_a_from_Jorgensen/BC_O2_Margin_updated.mat' , 'BC_O2_Margin_updated')
                txt_file = 'O2_muM_WOA2018_';
                txt_title = 'O_2 at the SWI (\muM)';
                levels = [0:10:310.0];
                limits = [0 300.0];
                
            case 3  % NO3
                BC_NO3_Margin_updated = BCond_BW_hr;
                save('BC_NO3_Margin_updated.mat' , 'BC_NO3_Margin_updated')
                save('./BC_calc_a_from_Jorgensen/BC_NO3_Margin_updated.mat' , 'BC_NO3_Margin_updated')
                txt_file = 'NO3_muM_WOA_';
                txt_title = 'NO_3 at the SWI (\muM)';
                levels = [0:2:52.0];
                limits = [0 50.0];
        end
        
        fig_SWI_BCond(bc) = figure;
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        [C,h] = m_contourf(long, lat,  BCond_BW_hr, levels);
        set(h,'LineColor','none')
        title(txt_title)
        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig_SWI_BCond(bc),'-depsc2', ['./BC_calc_a_from_Jorgensen/SWI_' txt_file, str_date '_HR_FILLED.eps']);
        
        % fig_SWI_BCond(bc) = figure;
        % set(gca,'FontSize',30)
        % pcolor(long, lat, BCond_BW_hr);
        % title(txt_title)
        % hold on;
        % shading interp;
        % %contour(long, lat, toc,'LineColor','k')
        % colorbar ('horizontal')
        % colormap(parula)
        % %caxis([-5000 0.0])
        % xlabel('Longitude')
        % ylabel('Latitude')
        % print(fig_SWI_BCond(bc),'-depsc2', ['./BC_calc_a_from_Jorgensen/SWI_' txt_file, str_date '_HR_FILLED.eps']);
        
    end
end
