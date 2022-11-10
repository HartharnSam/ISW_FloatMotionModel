% Plot for a DJL wave, the current flow speed and front and rear of float
% local flow speed

clc; clearvars; close all;
isRecalculateDJL = false;

LfLambdaList = 1:0.1:1.5; %0.03 0.5 1
%LfLambdaList = [0.03 0.3 1 .2163];
LfLambdaList = [0.01 0.03 0.1 0.3 0.5 1 1.2 1.5 2 3];
LfLambdaList = [.1523];

for j = 1:length(LfLambdaList)

LfLambdaList = [.1523];
    LfLambda = LfLambdaList(j);
    if isRecalculateDJL
        L = 14; H = 0.3;
        A = 1e-4;
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
        NXlist = [64    128 256 512 2048];
        NZlist = [32    64  128 256 1024];
        ddlist = [0.01  0.01 0.01 0.01  0.01];
        epsilonlist =[1e-4 1e-4 1e-4 1e-4 1e-5];
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
        clearvars -except c uwave x z density L wavelength max_u_par LfLambda j isRecalculateDJL
        save('../../02_Raw_data/DJL_Wave1', 'x', 'uwave', 'c', 'z', 'density', 'L', 'wavelength');
    else
        clearvars -except c uwave x z density L wavelength max_u_par LfLambda j isRecalculateDJL
        %load('../../02_Raw_data/DJL_Wave1');
        load('DJL')
        c = DJL.WaveC;
        uwave = DJL.u;

    end
    %%
    % Set up timestepping
    t1 = 0; t2 = 100;
    dt = 1/10;
    t = t1:dt:t2;
    % Set up a moving frame of reference for the DJL solution, set the starting
    % wave location as x=0
    x_cur = x' - c*t + L/2;

    % Calculate a moving frame of reference u profiles
    for ii = 1:length(t)
        u(:, ii) = interp1(x, uwave(end, :), x_cur(:, ii), 'linear', 'extrap');
    end
    partial_u = u;%(:, 1:1/dt:end);
    partial_t = t;%(1:1/dt:end);

    figure(1)
    pcolor(x, partial_t, partial_u'); %and plot
    cmocean('balance', 'pivot', 0);
    clrbar = colorbar;

    %% Parse and run model
    Flow.U_flow = u;
    Flow.timestep = t(2)-t(1);
    Flow.x = x';
    Flow.rho_0 = 1029;

    Particle.r = LfLambda*wavelength/2;

    Particle.StartLoc = wavelength + Particle.r +.5; % Start the particle just outside the wave's reach
    Particle.C_d = 170;
    Particle.rho_f = 910;
    Particle.Shape = 'Rectangle';

    [particle, fluid_u] = FloatMotionModel(Flow, Particle, 'basic');

    %% Plot
    figure(1);
    hold on
    plot(particle.x, t);

    % Figure 2
    figure(2);
    subplot(3, 1, 1)
    plot(t, particle.x);
    ylabel('x')
    subplot(3, 1, 2)
    plot(t, particle.u);
    ylabel('u')
    try
        subplot(3, 1, 3);
        plot(t, particle.dudt);
        ylabel('du_{}dt')
    end
    if mod(j, 1) == 0
        close all;
        %Figure 3
        figure(3);
        tiledlayout(2, 1)
        nexttile
        plot(t, particle.u/c, 'k-');
        hold on
        for i = 1:length(t)
            front_ind = nearest_index(x, particle.x(i)+Particle.r);
            rear_ind = nearest_index(x, particle.x(i)-Particle.r);
            front_u(i) = u(front_ind, i);
            rear_u(i) = u(rear_ind, i);
        end
        plot(t, rear_u/c, '-r');
        plot(t, front_u/c, 'b');
        max_u = .5;
        yline(0, '-','Color', [1 1 1]*.3);
        ylim([-max_u max_u])
        xlim([0 35])
        ylabel('$u/c_{isw}$', 'interpreter', 'latex')
        xticklabels([])
        legend('Float', 'Front Fluid', 'Rear Fluid', 'Location', 'best');
        title(['$L_f/\lambda = $ ', num2str(LfLambda)], 'interpreter', 'latex')

        % Add on difference in velocity part
        nexttile;
        plot(t, (particle.u-rear_u')/c, '-r');
        hold on
        plot(t, (particle.u - front_u')/c, '-b');
        ylim([-max_u max_u]);
        xlim([0 35])
        yline(0, '-','Color', [1 1 1]*.3);
        ylabel('$u_f - u(x) / c_{isw}$', 'interpreter', 'latex')
        xlabel('t (s)')

        figure_print_format(gcf, 15)
        fig = gcf;
        fig.Position = [681 331 560 467];
%        exportgraphics(gcf,['../../04_Output/06_SurfaceFlow/BasicFlowFloatModel_', num2str(LfLambda), '.eps'], 'ContentType', 'vector')
%        exportgraphics(gcf,['../../04_Output/06_SurfaceFlow/BasicFlowFloatModel_', num2str(LfLambda), '.png'])
        %exportgraphics(gcf,['../../04_Output/06_SurfaceFlow/FloatModels/BasicFlowFloatModel_', num2str(LfLambda), '.eps'], 'ContentType', 'vector')
        exportgraphics(gcf,['BasicFlowFloatModel_', num2str(LfLambda), '.png'])
    end
    

    max_u_par(j) = (max(particle.u)/c);

end
