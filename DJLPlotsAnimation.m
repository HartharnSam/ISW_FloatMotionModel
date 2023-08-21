% DJLPlotsAnimation - Produces an animation of the DJL Speed Plots function

clc; clearvars; close all;
LfLambda = 2.4/1.4712;
isRecalculateDJL = false;

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
    load('DJL');
end
%%
t1 = 0; t2 = 50;
dt = 1;
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

%% Parse and run model
Flow.u_flow = u;
Flow.timestep = t(2)-t(1);
Flow.x = x';
Flow.rho_1 = 1029;

Particle.r = LfLambda*wavelength/2;

Particle.StartLoc = wavelength + Particle.r; % Start the particle just outside the wave's reach
Particle.C_d = 1700;
Particle.rho_f = 910;
Particle.Shape = 'Rectangle';

[particle] = FloatMotionModel(Flow, Particle, 'basic');
tl = tiledlayout(3, 1, 'TileSpacing', 'tight');
ax1 = nexttile;
ax2 = nexttile;
ax3 = nexttile;
pDepth = -0.02;
fig = gcf;
fig.Units = 'centimeters';
fig.Position = [0 0 16 11.5];
%vid = VideoWriter('../../04_Output/06_SurfaceFlow/DJL_Animations.mp4', 'MPEG-4');
vid = VideoWriter('DJL_Animations.mp4', 'MPEG-4');

vid.FrameRate = 1;
open(vid)

for ii = 1:length(t)
    % Plot the xz plane
    axes(ax1)
    ut = interp2(x,z,  uwave, x_cur(:, ii), z, 'linear', 0);
    pcolor(x, z, ut/c);
    caxis(0.5*[-1 1]);
    cmocean('balance'); cb = colorbar; ylabel(cb, '$u/c$','interpreter', 'latex')
    xlim([0 7])
    hold on
    area([particle.x(ii)-Particle.r particle.x(ii)+Particle.r], pDepth*[1 1], 'FaceColor', 'w','EdgeColor', 'k');
    plot(particle.x(ii)-Particle.r, pDepth, 'or')
    plot(particle.x(ii)+Particle.r, pDepth, 'ob')
    plot(particle.x(ii), pDepth/2, 'xk')
    hold off
    ylabel('$z (m)$', 'interpreter', 'latex')
    yticks([-.3 0])
    ylim([-.3 0])
    % Plot the evolution of u's
    axes(ax2)
    plot(t(1:ii), particle.u(1:ii)/c, 'k-');
    hold on
    front_ind = nearest_index(x, particle.x(ii)+Particle.r);
    rear_ind = nearest_index(x, particle.x(ii)-Particle.r);
    front_u(ii) = u(front_ind, ii);
    rear_u(ii) = u(rear_ind, ii);

    plot(t(1:ii), rear_u/c, '-r');
    plot(t(1:ii), front_u/c, 'b');

    yline(0, '-', 'color', [.5 .5 .5])
    max_u = .6;
    ylim([-max_u max_u])
    xlim([t1 t2])
    hold off
    ylabel('$u/c_{isw}$', 'interpreter', 'latex')
    legend('Float', 'Fluid A', 'Fluid B','', 'Location', 'eastoutside');
    xticklabels([]);

    axes(ax3)
    plot(t(1:ii), (particle.u(1:ii)-rear_u')/c, '-r');
    hold on
    plot(t(1:ii), (particle.u(1:ii) - front_u')/c, '-b');
    yline(0, '-', 'color', [.5 .5 .5])
    max_u = .5;
    ylim([-max_u max_u])
    xlim([t1 t2])
    hold off
    ylabel('$(u_f-u(x))/c_{isw}$', 'interpreter', 'latex')
    xlabel('$t (s)$');

    drawnow;
    hold off
    pause(.1)
    figure_print_format(gcf, 18);
    dark_figure(gcf, [23 23 23])
    F = getframe(gcf);
    writeVideo(vid, F);
    completion(ii-t1, t2-t1);
end
close(vid);
