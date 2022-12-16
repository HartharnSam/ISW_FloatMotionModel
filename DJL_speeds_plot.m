%DJL_SPEEDS_PLOT - Plots the fluid speeds at point A and point B for a
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
isRecalculateDJL = false;
makePlots = true;
%LfLambdaList = logspace(-1.5229,0.7782, 20)';
%LfLambdaList = [0.03 0.3 3];
lambdas = logspace(-5, -3, 10);
lambdas = 1.4712;
LfLambdaList = 1.2/lambdas;
max_u_par = LfLambdaList*lambdas*NaN;
list_wavelength = lambdas*NaN;
list_amp = lambdas*NaN;
list_c = lambdas*NaN;

for i_lambda = 1:length(lambdas)

    for j = 1:length(LfLambdaList)
        LfLambda = LfLambdaList(j);
        if isRecalculateDJL && j == 1
            L = 14; H = 0.3;
            A = lambdas(i_lambda);
            verbose = 0; relax = 0.15;

            % Specify the general density profile which takes d_d as a second parameter
            a_d = 0.019/2;
            z0_d = 0.07;

            frho = @(z, d_d) 1-a_d*tanh((z+z0_d)/d_d);
            frhoz = @(z, d_d) -(a_d/d_d)*sech((z+z0_d)/d_d).^2;

            % The velocity profile (zero for this case) (m/s)
            Ubg = @(z) 0*z; Ubgz = @(z) 0*z; Ubgzz = @(z) 0*z;

            % Find solution
            start_time = clock;

            % Specify resolution and pycnocline parameter d_d
            NXlist = [64    128 256 512 1024 2048];
            NZlist = [32    64  128 256 256 512];
            ddlist = [0.01  0.01 0.01 0.01  0.01 .01];
            epsilonlist =[1e-4 1e-4 1e-4 1e-4 1e-4 1e-5];
            for ddindex = 1:length(ddlist)
                NX = NXlist(ddindex);
                NZ = NZlist(ddindex);

                d_d = ddlist(ddindex);

                rho = @(z) frho(z, d_d);
                rhoz = @(z) frhoz(z, d_d);
                epsilon = epsilonlist(ddindex);
                djles_refine_solution;
            end
            end_time = clock;

            djles_diagnostics;
            clearvars -except c uwave x z density L wavelength max_u_par LfLambda* j isRecalculateDJL  max_u_par makePlots lambdas wave_ampl list_* i_lambda
            save('../../02_Raw_data/DJL_Wave_tmp', 'x', 'uwave', 'c', 'z', 'density', 'L', 'wavelength', 'wave_ampl');
        else
            clearvars -except c uwave x z density L wavelength max_u_par LfLambda* j isRecalculateDJL max_u_par makePlots lambdas list_* i_lambda
            %load('../../02_Raw_data/DJL_Wave_tmp');
            load('DJL')
            c = DJL.WaveC;
            uwave = DJL.u;
            wave_ampl = -DJL.WaveAmp;
            %c = .107;
        end
        %%
        % Set up timestepping
        t1 = 0; t2 = 100;
        dt = 1/10;
        tim = t1:dt:t2;
        % Set up a moving frame of reference for the DJL solution, set the starting
        % wave location as x=0
        x_cur = x' - c*tim + L/2;

        % Calculate a moving frame of reference u profiles
        u = x_cur*NaN;
        for ii = 1:length(tim)
            u(:, ii) = interp1(x, uwave(end, :), x_cur(:, ii), 'linear', 'extrap');
        end
        partial_u = u;%(:, 1:1/dt:end);
        partial_t = tim;%(1:1/dt:end);

        %figure(1)
        %pcolor(x, partial_t, partial_u'); %and plot
        %cmocean('balance', 'pivot', 0);
        %clrbar = colorbar;

        %% Parse and run model
        Flow.u_flow = u;
        Flow.timestep = tim(2)-tim(1);
        Flow.x = x';
        Flow.rho_1 = 1029;

        Particle.r = LfLambda*wavelength/2;

        Particle.StartLoc = wavelength + Particle.r +.5; % Start the particle just outside the wave's reach
        Particle.C_d = 170;
        Particle.rho_f = 910;
        Particle.Shape = 'Rectangle';

        [particle, fluid_u] = FloatMotionModel(Flow, Particle, 'basic');

        %% Plot
        if makePlots
            figure(1);
            hold on
            plot(particle.x, tim);

            % Figure 2
            figure(2);
            subplot(3, 1, 1)
            plot(tim, particle.x);
            ylabel('x')
            subplot(3, 1, 2)
            plot(tim, particle.u);
            ylabel('u')
            try
                subplot(3, 1, 3);
                plot(tim, particle.dudt);
                ylabel('du_{}dt')
            end
        end
        if mod(j, 1) == 0
            if makePlots
                close all;
                %Figure 3
                figure(3);
                tiledlayout(3, 1)
                nexttile
                plot(tim, particle.u/c, 'k-');
                hold on
            end
            front_u = tim*NaN; rear_u = tim*NaN;
            for i = 1:length(tim)
                front_ind = nearest_index(x, particle.x(i)+Particle.r);
                rear_ind = nearest_index(x, particle.x(i)-Particle.r);
                front_u(i) = u(front_ind, i);
                rear_u(i) = u(rear_ind, i);
            end
            if makePlots
                plot(tim, rear_u/c, '-r');
                plot(tim, front_u/c, 'b');
                max_u = .6;
                yline(0, '-','Color', [1 1 1]*.3);
                ylim([-max_u max_u])
                xlim([0 70])
                ylabel('$u/c_{isw}$', 'interpreter', 'latex')
                xticklabels([])
                title(['$L_f/\lambda = $ ', num2str(LfLambda)], 'interpreter', 'latex')
                hold on
                %xline([6 11 20.5 31 37])
                %legend('Float', 'Fluid A', 'Fluid B','', '', '', '', '', 'Location', 'eastoutside');


                % Add on difference in velocity part
                nexttile;
                plot(tim, (particle.u-rear_u')/c, '-r');
                hold on
                plot(tim, (particle.u - front_u')/c, '-b');
                ylim([-max_u max_u]);
                xlim([0 70])
                yline(0, '-','Color', [1 1 1]*.3);
                ylabel('$u_f - u(x) / c_{isw}$', 'interpreter', 'latex')
                xticklabels([])

                nexttile
                plot(tim, particle.x-particle.x(1), '-r');
                xlim([0 70]); ylim([0 1]);
                yline(0, '-', 'Color', [1 1 1]*.3);
                ylabel('$x_f$','interpreter', 'latex');
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
        end


        max_u_par(j, i_lambda) = (max(particle.u)/c);
        list_wavelength(i_lambda) = wavelength;
        list_amp(i_lambda) = -wave_ampl;
        list_c(i_lambda) = c;
        
    end
end

max_uf_c.data = max_u_par;
max_uf_c.lambda = list_wavelength;
max_uf_c.LfLambda = LfLambdaList;
max_uf_c.amp = list_amp;
max_uf_c.c = list_c;
max_uf_c.APE = lambdas;
save('max_uf_c.mat', 'max_uf_c')