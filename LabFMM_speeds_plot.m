%LAB_PTM_SPEEDS_PLOT - Plots the fluid speeds at point A and point B for a
%float moving with the fluid according to the FloatMotionModel
%
% Other m-files required: FloatMotionModel, cmocean, djles,
% figure_print_format, nearest_index
% Subfunctions: none
% MAT-files required: none
%
% See also:
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 10-Sep-2022; Last revision: 13-Dec-2022
% MATLAB Version: 9.10.0.1602886 (R2021a)

digiflowstartup
clc; clearvars; close all;
makePlots = true;
% Need to run for 090322, 141122, 141222
c = 0.107;
%LfLambda = %
        %%
        % Set up timestepping
filename = {'./CamA/piv_ts.dfi', './CamB/piv_ts.dfi'};
isFillMissing = false; 
t_start = 0; 

if size(filename, 2)==1 %% Just be nice and load in one set of data
    im = dfireadvel(filename{1});
    u = im.cdata(:, :, 1);
    u(u == 0) = NaN; % Remove anomalous numbers
    u = flip(u, 1)';
    x = im.x;
    x = x(1, :)';
    times = flip(im.y(:, 1));

elseif size(filename, 2) == 2 %% Merge two timeseries
    im1 = dfireadvel(filename{1});
    im2 = dfireadvel(filename{2});
    grid_1 = dfi_grid_read(im1);
    grid_2 = dfi_grid_read(im2);

    cutoff = 8; % Amount to cut off the edges, which appears to be an effect on PIV images
    im1.cdata(:, [1:cutoff end-cutoff:end], :) = NaN;
    im2.cdata(:, [1:cutoff end-cutoff:end], :) = NaN;
    
    xmin= min(min(grid_1.x), min(grid_2.x));
    xmax= max(max(grid_1.x), max(grid_2.x));
    tmin= min(min(grid_1.y), min(grid_2.y));
    tmax= max(max(grid_1.y), max(grid_2.y));
    
    new_x = xmax:grid_1.dx:xmin;
    new_t = tmin:grid_1.dy:tmax;
    [newX, newT] = meshgrid(new_x, new_t);

    new_1_data = interp2(grid_1.X, grid_1.Y, im1.cdata(:, :, 1), newX, newT);
    new_2_data = interp2(grid_2.X, grid_2.Y, im2.cdata(:, :, 1), newX, newT);
    u = new_1_data;
    u(~isnan(new_2_data)) = new_2_data(~isnan(new_2_data));
    u = u';
    times = new_t';
    x = new_x';
else
    error('Length of filename wrong')
end

if isFillMissing
    u = fillmissing(u, 'linear', 'EndValues', 'none')'; %Re-fill those values
end

%% Cut down to the requested timings
time_index = nearest_index(times, t_start):length(times);
times = times(time_index);
u = u(:, time_index);

Flow.U_flow = u;
Flow.timestep = times(2)-times(1);
Flow.x = x;

        %figure(1)
        %pcolor(x, partial_t, partial_u'); %and plot
        %cmocean('balance', 'pivot', 0);
        %clrbar = colorbar;

        %% Parse and run model
        Flow.u_flow = u;
        Flow.timestep = times(2)-times(1);
        Flow.x = x;
        Flow.rho_1 = 1029;

        Particle.r = .1/2;

        Particle.StartLoc = 4.2201; % Start the particle just outside the wave's reach
        Particle.C_d = 170;
        Particle.rho_f = 910;
        Particle.Shape = 'Circle';

        [particle, fluid_u] = FloatMotionModel(Flow, Particle, 'basic');

        %% Plot
        if makePlots
            figure(1);
            hold on
            plot(particle.x, times);

            % Figure 2
            figure(2);
            subplot(3, 1, 1)
            plot(times, particle.x);
            ylabel('x')
            subplot(3, 1, 2)
            plot(times, particle.u);
            ylabel('u')
            try
                subplot(3, 1, 3);
                plot(times, particle.dudt);
                ylabel('du_{}dt')
            end
        end
            if makePlots
                close all;
                %Figure 3
                figure(3);
                tiledlayout(2, 1)
                nexttile
                plot(times, particle.u/c, 'k-');
                hold on
            end
            front_u = times*NaN; rear_u = times*NaN;
            for i = 1:length(times)
                front_ind = nearest_index(x, particle.x(i)+Particle.r);
                rear_ind = nearest_index(x, particle.x(i)-Particle.r);
                front_u(i) = u(front_ind, i);
                rear_u(i) = u(rear_ind, i);
            end
            if makePlots
                plot(times, rear_u/c, '-r');
                plot(times, front_u/c, 'b');
                max_u = .6;
                yline(0, '-','Color', [1 1 1]*.3);
                ylim([-max_u max_u])
                xlim([0 70])
                ylabel('$u/c_{isw}$', 'interpreter', 'latex')
                xticklabels([])
%                title(['$L_f/\lambda = $ ', num2str(LfLambda)], 'interpreter', 'latex')
                hold on
                %xline([6 11 20.5 31 37])
                %legend('Float', 'Fluid A', 'Fluid B','', '', '', '', '', 'Location', 'eastoutside');


                % Add on difference in velocity part
                nexttile;
                plot(times, (particle.u-rear_u)/c, '-r');
                hold on
                plot(times, (particle.u - front_u)/c, '-b');
                ylim([-max_u max_u]);
                xlim([0 70])
                yline(0, '-','Color', [1 1 1]*.3);
                ylabel('$u_f - u(x) / c_{isw}$', 'interpreter', 'latex')
                %xticklabels([])
% 
%                 nexttile
%                 plot(times, particle.x-particle.x(1), '-r');
%                 xlim([0 70]); ylim([0 1]);
%                 yline(0, '-', 'Color', [1 1 1]*.3);
%                 ylabel('$x_f$','interpreter', 'latex');
 xlabel('t (s)')

                figure_print_format(gcf, 18)
                fig = gcf;
                fig.Units = 'centimeters';
                fig.Position = [0 0 14 10.5];
                exportgraphics(gcf, 'DJL_FluidSpeeds.png')
            end
            %        exportgraphics(gcf,['../../04_Output/06_SurfaceFlow/BasicFlowFloatModel_', num2str(LfLambda), '.eps'], 'ContentType', 'vector')
            %        exportgraphics(gcf,['../../04_Output/06_SurfaceFlow/BasicFlowFloatModel_', num2str(LfLambda), '.png'])
            %exportgraphics(gcf,['../../04_Output/06_SurfaceFlow/FloatModels/BasicFlowFloatModel_', num2str(LfLambda), '.eps'], 'ContentType', 'vector')
            %dark_figure(gcf, [23 23 23])
            %export_fig(gcf,['BasicFlowFloatModel_', num2str(LfLambda), '.png'], '-dpng')

            %  exportgraphics(gcf,['BasicFlowFloatModel_', num2str(LfLambda), '.png'])


        max_u_par(j, i_lambda) = (max(particle.u)/c);
        list_wavelength(i_lambda) = wavelength;
        list_amp(i_lambda) = -wave_ampl;
        list_c(i_lambda) = c;
        

max_uf_c.data = max_u_par;
max_uf_c.lambda = list_wavelength;
max_uf_c.LfLambda = LfLambdaList;
max_uf_c.amp = list_amp;
max_uf_c.c = list_c;
max_uf_c.APE = lambdas;
save('max_uf_c.mat', 'max_uf_c')