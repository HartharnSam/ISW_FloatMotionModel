function [particle_t, particle_x, particle_u, u_data] = ParticleTrackModel(C_d, rho_f,x_start_loc, isPlot)
%ADVANCED_PTM - Model of particle motion using SPINS input where:
%   Force = 1/2 rho_0 U^2 A C_d - where C_d is drag coefficient, U is
%   relative velocity, A is object area
%   du_dt = F/m , so A/m = 1/rho_f
% Inputs:
%    C_d - Test drag coefficient
%    rho_f - Density of float
%    x_start_loc - Initial location of float (index of x)
%    isPlot - Boolean to plot at end
%
% Outputs:
%    particle_t - Vector of times
%    particle_x - Vector of the particle's x location at each corresponding
%    time
%    particle_u - Vector of the particle's velocity (u) at each
%    corresponding time
%    u_data -
% TODO: Integrate into the advanced_PTM model o keep both plotting &
% calculation tools together
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 19-Jan-2022; Last revision: 19-Jan-2022
% MATLAB Version: 9.10.0.1602886 (R2021a)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

%clc; clearvars; close all;
spinsstartup;

orig_dir = cd;
addpath(orig_dir);
ICW_Path = ['C:\Users\', getenv('USERNAME'), '\OneDrive - Newcastle University\02_PhD_Project\04_Ice_Covered_Waters\'];
cd([ICW_Path, '\02_Raw_data\SurfaceFlowSPINS\03_170421']);
%cd(['C:\Users\', getenv('USERNAME'), '\OneDrive - Newcastle University\02_PhD_Project\03_ShoalingStratification\02_Raw_data\Model\NovakTank\081020_44']);

%% Set Parameters
% For changing Pycnocline Location
filename_ends = { '../03_170421'};%,'../03_170421', '../08_220421'};
%file_text = { '$pyc_{loc} = 0.05$', '$pyc_{loc} = 0.07*$','$pyc_{loc} = 0.09$'};
%save_fnm = 'pyclocVarying'

% For a slope
%filename_ends = { '../250720_31'};%,'../03_170421', '../08_220421'};

savefnm = [orig_dir, 'BAMC.mp4'];
z_coord = -.01;
x_start_loc = [2000 2500 3000]/2;
mode = 'animation';
background = 'hovmoller';
clims = .1*[-1 1];
%
% filename_ends = {'../05_190421','../03_170421',  '../06_200421'};
% file_text = {'$$\Delta\rho = 0.009$$', '$$\Delta \rho = 0.019*$$','$$\Delta\rho = 0.028$$'};
% save_fnm = 'delrhoVarying';
%% Start Script
fig = figure;
cd(filename_ends{1});
[~, x, data] = SPINS_hovmoller_z('u',z_coord, []);
close(fig);
t = 1:size(data, 2);
t = 1:60;
x_pos = nan(length(x_start_loc), length(t)+1);
u_t = x_pos;
x_pos(:, 1) = x(x_start_loc);

%% New x locations

if strcmpi(mode, 'animation')
    % Get grid ready for animation panel
    gd.x = xgrid_reader;
    x_inds = find(gd.x(:, 1)>3 & gd.x(:, 1)<7);
    gd.x = gd.x(x_inds, :);
    gd.z = zgrid_reader;
    gd.z = gd.z(x_inds, :);
    % And the video ready
    vid = VideoWriter(savefnm, 'MPEG-4');
    vid.FrameRate = 1;
    open(vid);
    n_panels = 2;
else
    n_panels = 1;
end
if strcmpi(background, 'hovmoller')
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
TypeParameter.Type  = 'model';
TypeParameter.Model = 'advanced';
isPlot = false;
[particle_t, x_pos] = advanced_PTM(TypeParameter, C_d, rho_f, x_start_loc)

for i = t
%     for j = 1:length(x_start_loc)
%         u_t(j, i) = interp1(x, data(:, i), x_pos(j, i));
%         x_pos(j, i+1) = x_pos(j, i)+u_t(j, i);
%     end
    subaxis(n_panels, 1, n_panels, 'MarginRight', .2)
    
    if strcmpi(background, 'hovmoller')
        pcolor(x, t([1 1:i]), data(:, [1 1:i])');
        shading flat
        caxis(clims)
        cmocean('balance')
        hold on
        try
            p1 = plot(wave_center(1:i), time(1:i), 'k:');
            plot(wave_center(1:i) + wavelength_right(1:i), time(1:i),'k--')
            p2 = plot(wave_center(1:i) - wavelength_left(1:i), time(1:i), 'k--');
        catch
            p1 = plot(wave_center, time, 'k:');
            plot(wave_center + wavelength_right, time,'k--')
            p2 = plot(wave_center - wavelength_left, time, 'k--');
            
        end
        ax1 = gca;
        ax1_pos = ax1.Position;
        c = colorbar;
        ylabel(c, 'u')
        ax1.Position = ax1_pos;
        %legend([p1 p2], 'Wave Center', '$\pm 1$ Wavelength')
        
    end
    p3 = plot((x_pos(:, 1:end-1))', [t; t; t]');
    ylabel('Time (s)')
    xlabel('X Position (m)')
    xlim([3 6.6])
    ylim([0 max(t)])
    hold off
    if strcmpi(mode, 'animation')
        ax2 = subaxis(n_panels, 1, 1, 'MarginRight', .2);
        temp_data = spins_reader_new('u', i-1, x_inds, []);
        pcolor(gd.x, gd.z, temp_data),shading flat
        caxis(clims)
        cmocean('balance');
        hold on
        plot(x_pos(:, i), z_coord, 'x');
        
        ax2_pos = ax2.Position;
        c = colorbar;
        ylabel(c, '$u / c$');
        xlim([3 6.6])
        hold off
        
        ax2.Position = ax2_pos;
        
    end
    figure_print_format(gcf);
    dark_figure(gcf, [23 23 23]);
    
    if strcmpi(mode, 'animation')
        F = getframe(gcf);
        writeVideo(vid, F);
    end
end
if ~strcmpi(mode, 'animation')
    print([savefnm, '.png'], '-dpng');
    print([savefnm, '.eps'], '-depsc');
    
else
    close(vid)
end

cd(orig_dir);