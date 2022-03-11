function plot_ATM(flow, particle, u, PlotType, plotchars, DataType)
% flow.x, flow.z, flow.t, u, particle
% if model or lab
% xlimits
% output_fname
xlimits = plotchars.xlimits;
output_fname = plotchars.output_fname;
if strcmp(DataType, 'model')
    model = true;
    lab = false;
elseif strcmpi(DataType, 'lab')
    lab = true;
    model = false;
end

switch PlotType
    case 'animation'
        vid = VideoWriter('output_fname', 'MPEG-4');
        vid.FrameRate = 1/(flow.t(2)-flow.t(1));
        open(vid);
        
        tiledlayout(2, 1);
        ax1 = nexttile;
        ax2 = nexttile;
        
        xinds = find(flow.x(:, 1) > xlimits(1) & flow.x(:, 1) < xlimits(2));
        clims = [-1 1].*abs(max(u(:)));
        clrorder = custom_color_orders('set_colormap', 'Art', 'Juarez', size(particle, 2));
        
        if model
            x = flow.x(xinds, :);
            z = flow.z(xinds, :);
        end
        
        if isfield('plotchars', 'wave_props') && plotchars.wave_props
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
        
        for ii = 1:length(flow.t)
            time = flow.t(ii);
            
            if model
                temp_data = spins_reader_new('u', time, xinds, []);
            elseif lab
                frm_fname = strrep('piv_####.dfi', '####', sprintf('%04d', time));
                tmpim = dfireadvel(frm_fname);
                temp_data = tmpim(xinds, :, 1);
            end
            
            axes(ax1);
            pcolor(ax1, flow.x(xinds, :), flow.x(xinds, :), temp_data);
            shading flat
            caxis(ax1, clims)
            cmocean('balance');
            hold(ax1, 'on');
            plot(ax1, particle.x(ii, :), -.05, 'o', 'MarkerFaceColor', 'auto', 'MarkerSize', 8);
            c = colorbar(ax1);
            ylabel(c, 'u (m/s)')
            xlim(xlimits)
            hold off;
            
            axes(ax2);
            pcolor(ax2, flow.x, flow.t([1 1:ii]), u(:, [1 1:ii])');
            hold(ax2, 'on');
            shading flat
            caxis(ax2, clims);
            cmocean('balance');
            
            if isfield('plotchars', 'wave_props') && plotchars.wave_props
                plot(ax2, wave_center(1:ii), time(1:ii), 'k:');
                plot(ax2, wave_center(1:ii) + wavelength_right(1:ii), time(1:ii),'k--')
                plot(ax2, wave_center(1:ii) - wavelength_left(1:ii), time(1:ii), 'k--');
            end
            c = colorbar(ax2);
            ylabel(c, 'u');
            
            custom_color_orders('set_colormap', 'Art', 'Juarez', length(x_start_ind));
            plot(ax2, particle.x(1:ii, :), flow.t(1:ii), '-');
            xlim(ax2, xlimits)
            ylim(ax2, [flow.t(1) flow.t(end)]);
            hold off;
            
            %dark_figure(gcf, [23 23 23]);
            figure_print_format(gcf);
            
            F = getframe(gcf);
            writeVideo(vid, F);
        end
        close(vid);
        
        %%
    case 'ParticleHovmoller'
        fig = figure;
        pcolor(flow.x, flow.t([1 1:ii]), u(:, [1 1:ii])');
        hold('on');
        shading flat
        caxis(clims);
        cmocean('balance');
        
        if wave_props
            plot( wave_center(1:ii), time(1:ii), 'k:');
            plot( wave_center(1:ii) + wavelength_right(1:ii), time(1:ii),'k--')
            plot( wave_center(1:ii) - wavelength_left(1:ii), time(1:ii), 'k--');
        end
        c = colorbar;
        ylabel(c, 'u');
        
        custom_color_orders('set_colormap', 'Art', 'Juarez', length(x_start_ind));
        plot( particle.x(1:ii, :), flow.t(1:ii), '-');
        xlim( xlimits)
        ylim([flow.t(1) flow.t(end)]);
        hold off;
        set(gca, 'XDir', 'reverse');
        
        %dark_figure(gcf, [23 23 23]);
        figure_print_format(gcf);
        exportgraphics(gcf, output_fname);
    otherwise
        error('Unknown plot type')
        
end



