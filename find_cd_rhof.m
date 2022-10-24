%fit_ISWPTM
clc; clearvars; close all; digiflowstartup;

m_path = mpath;
%cd(m_path);
addpath(m_path);

%% Get the lab Particle Data
load('./CamC/ptv_tracks.mat', 'ptv'); % Load in the data
im = dfireadvel('CamC/output_0000.dfi'); % Load in a frame to get WCS adjustment
Grid = dfi_grid_read(im);

%TODO Not just the longest, but one that starts earlier
for i = 1:ptv.n_particles
    if ~isempty(ptv.data{i})
        locations_tmp = interp1([1 Grid.nx], (Grid.x), ptv.data{i}(:, 32));
        locations_tmp = locations_tmp(~isnan(locations_tmp)); % Remove nans
    end
    if exist('locations_tmp', 'var')
        break
    end
end

times = (0:ptv.n_timesteps-1)+1; % Time indices
times = times(~isnan(locations_tmp))/30; % Remove nans
locations = locations_tmp(~isnan(locations_tmp))'; % Remove nans

%% Make first guess at least squares fit
coefs_0 = 10;%, 910, locations(1)]; % [C_d, rho_f, x_start];
t_start = nearest_index(times, 3.1667); 
t_end = nearest_index(times, 21.433);
times = times(t_start:t_end);
locations = locations(t_start:t_end);
disp(num2str(locations(1)))
coefs_1 = lsqcurvefit(@test_FMM, coefs_0, times, locations);

[locations_fitted] = test_FMM(coefs_1, times);
fprintf('\n C_d =  %4.2f \n', coefs_1(1));

%% Plot the data
figure
plot(locations, times, 'k-');  % Plot the lab data
hold on
plot(locations_fitted, times, 'b-')

% Make plot nice
xlabel('x (m)');
ylabel('t (s)');
l = legend('Lab Measured', ['Fit : $C_d$:', num2str(coefs_1(1))]);
l.Interpreter = 'latex';
ax = gca; ax.XDir = 'reverse';
figure_print_format(gcf, 18);

fig = gcf; 
fig.Color = [23 23 23]/255;
ax.XColor = 'w';
ax.YColor = 'w';
c.YColor = 'w';

fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 18 10];
wysiwyg;

%export_fig(fig, '../../../06_Communication/BAMC/LabFittedParticles.png', '-dpng');


function x = test_FMM(coefs_0, times)
filename = {'./CamA/piv_ts.dfi', './CamB/piv_ts.dfi'};
t_start = times(1); 
t_end = times(end);
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

time_index = nearest_index(times, t_start):nearest_index(times, t_end);
times = times(time_index);
u = u(:, time_index);

Flow.U_flow = u;
Flow.timestep = times(2)-times(1);
Flow.x = x;

Particle.StartLoc = 4.7925;
Particle.C_d = coefs_0;
Particle.rho_f = 910;

[particle, ~] = FloatMotionModel(Flow, Particle, 'advanced');
x = particle.x';

end
