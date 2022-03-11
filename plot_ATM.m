if isPlot && strcmpi(TypeParameter.PlotType, 'animation')
        if ii == 1
            [xgrid, zgrid] = spinsgrid2d;
            vid = VideoWriter('aPTM.mp4', 'MPEG-4');
            vid.FrameRate = 1;
            open(vid);
            tlayout = tiledlayout(2, 1);
            ax1 = nexttile;
            ax2 = nexttile;
            xinds = find(xgrid(:, 1) > xlimits(1) & xgrid(:, 1) < xlimits(2));
            clims = [-1 1].*abs(max(u(:)));
            
            try
                load wave_characteristics.mat wave_center wavelength_* time
            catch
                characterize_wave;
                load wave_characteristics.mat wave_center wavelength_* time
            end
            wave_center(wave_center == 0) = nan;
            wavelength_left(wavelength_left == 0) = nan;
            wavelength_right(wavelength_right == 0) = nan;
        end
        
        temp_data = spins_reader_new('u', ii-1, xinds, []);
        axes(ax1);
        
        pcolor(ax1, xgrid(xinds, :), zgrid(xinds, :), temp_data);
        shading flat
        caxis(ax1, clims)
        cmocean('balance');
        clrorder = custom_color_orders('set_colormap', 'Art', 'Juarez', length(x_start_ind));
        hold on
        plot(particle_x(ii, :), z_coord, 'o', 'MarkerFaceColor', 'auto', 'MarkerSize', 8);
        c = colorbar;
        ylabel(c, 'u')
        xlim(xlimits)
        hold off;
        
        axes(ax2);
        pcolor(X, times([1 1:ii]), u(:, [1 1:ii])');
        hold on
        shading flat
        caxis(clims)
        cmocean('balance')
        plot(wave_center(1:ii), time(1:ii), 'k:');
        plot(wave_center(1:ii) + wavelength_right(1:ii), time(1:ii),'k--')
        plot(wave_center(1:ii) - wavelength_left(1:ii), time(1:ii), 'k--');
        c = colorbar;
        custom_color_orders('set_colormap', 'Art', 'Juarez', length(x_start_ind));
        plot(particle_x(1:ii, :), times(1:ii), '-');
        
        ylabel(c, 'u')
        xlim(xlimits)
        ylim([0 length(times)-1]);
        hold off;
        
        dark_figure(gcf, [23 23 23]);
        figure_print_format(gcf);
        
        F = getframe(gcf);
        writeVideo(vid, F);
end
    

if isPlot && strcmpi(TypeParameter.PlotType, 'animation')
    close(vid);
end
if isPlot
    figure
    pcolor(X, time', u'); shading flat;
    caxis([-.1 .1]);
    newbluewhitered;
    %axis tight
    
    set(gca, 'XDir', 'reverse');
    
    %     plot(particle_x, particle_t, 'b-');
    %     plot(basic_particle_x, particle_t, 'r-');
    %
    %
    %     figure
    %     subaxis(2, 1, 1)
    %     plot(particle_t, particle_x, 'k-');
    %     hold on
    %     plot(particle_t, basic_particle_x, 'b-');
    %     subaxis(2, 1, 2)
    %     plot(particle_t, particle_u);
    %     hold on
    %     plot(particle_t, u_t);
end


pcolor(x, [0:65], u(:, 1:66)'); shading flat; caxis([-.1 .1]); newbluewhitered; cbar = colorbar
hold on;
plot(particle_x, particle_t, 'k-')
xlim([4 6]); ylim([30 60])

exportgraphics(gcf, 'D:\OneDrive - Newcastle University\02_PhD_Project\04_Ice_Covered_Waters\06_Communication\BAMC\ParticleTracks.png')

ax = gca;
ax.Position = [0 0 1 1];
ax.XColor = 'none';
ax.YColor = 'none';
delete cbar

fig = gcf; 
fig.Position(4) = fig.Position(3);
exportgraphics(gcf, 'D:\OneDrive - Newcastle University\02_PhD_Project\04_Ice_Covered_Waters\06_Communication\BAMC\ParticleTracks_thumb.png')



%%
%%
[particle_t, particle_x, particle_u, u, x] = advanced_PTM(TypeParameter, C_d, rho_f, x_start_loc, isPlot, u_0);
pcolor(particle_t, x, u); shading flat; caxis([-.1 .1]);
newbluewhitered;
% colorbar;
hold on
plot(particle_t, particle_x);

plot_ptv_tracks('x', '030222');