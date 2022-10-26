%FIT_CD
% Fit C_d
clc; clearvars; close all; digiflowstartup;

m_path = mpath;
addpath(m_path);

% First guess at times to start and end
piv = dfireadvel('CamA/piv_ts.dfi');
t_start = piv.yOriginWorld;
t_end = piv.yOriginWorld + ((piv.ny-1)*piv.yWorldPerPixel);

%% Get the lab Particle Data
load('./CamC/ptv_tracks_compiled.mat', 'ptv'); % Load in the data
im = dfireadvel('CamC/output_0000.dfi'); % Load in a frame to get WCS adjustment
Grid = dfi_grid_read(im);
times = (0:ptv.n_timesteps-1)+1; % Time indices
fin = 1;

particles = find(~cellfun(@isempty, ptv.data));
locations = cell(1, length(particles));
times_all = locations;
t_starts = NaN(1, length(particles));
t_ends = t_starts; loc_starts = t_starts;
long_tracks = false(1, length(particles));

for i = 1:length(particles)
    j = particles(i);
    locations_tmp = interp1([1 Grid.nx], (Grid.x), ptv.data{j}(:, 32));
    locations{i} = locations_tmp(~isnan(locations_tmp))'; % Remove nans
    times_all{i} = times(~isnan(locations_tmp))/30;

    %% Make first guess at least squares fit
    coefs_0 = 10;
    t_starts(i) = nearest_index(times_all{i}, t_start);
    t_ends(i) = nearest_index(times_all{i}, t_end);
    times_all{i} = times_all{i}(t_starts(i):t_ends(i));
    locations{i} = locations{i}(t_starts(i):t_ends(i));
    loc_starts(i) = locations{i}(1);
    long_tracks(i) = length(locations{i})>1;

    figure(fin);
    hold on
    plot(times_all{i}, locations{i}, 'k-')
end
% Cut down tracks to length>1 tracks only
times_all = times_all(long_tracks);
locations = locations(long_tracks);
t_starts = t_starts(long_tracks);
t_ends = t_ends(long_tracks);
loc_starts = loc_starts(long_tracks);
fin = fin+1;

%% Do Least Squares Fit to find drag coefficient
function_args.tstarts = t_starts;
function_args.t_ends = t_ends;
function_args.StartLoc = loc_starts;
[coefs_1, R2,~,~,~,~,~, argsIn] = lsqcurvefit_ul(@test_FMM, coefs_0, times_all, locations, function_args);

%%
num_times = 0;
for ii = 1:length(times_all)
    num_times = num_times+length(times_all{ii});
end
xdata = NaN(num_times, 1); ydata = xdata;
end_index = zeros(1, length(times_all)+1);

for ii = 1:length(times_all)
    indexes = end_index(ii)+(1:length(times_all{ii}));
    end_index(ii+1) = indexes(end);
    xdata(indexes) = times_all{ii};
    ydata(indexes) = locations{ii};
end
times_2 = xdata;
locations_2 = ydata;
%[locations_fitted] = test_FMM(coefs_1, times_2, argsIn);
fprintf('\n C_d =  %4.2f \n', coefs_1(1));

%% Test the residuals for various coefs_1
coefs_0 = [logspace(-2, 4,100)];
R2 = coefs_0*NaN;

parfor ii = 1:length(coefs_0)
    locations_tmp = test_FMM(coefs_0(ii), times_2, argsIn);
    R2(ii) = sum((locations_2-locations_tmp).^2);
end
clf;
plot(coefs_0, R2, 'xk');

%% Plot the data
figure
[locations_fitted] = test_FMM(coefs_1(end), times_2, argsIn);

plot(locations_2, times_2, 'kx');  % Plot the lab data
hold on
plot(locations_fitted, times_2, 'b.')

% Make plot nice
xlabel('x (m)');
ylabel('t (s)');
l = legend('Lab Measured', ['Fit : $C_d$:', num2str(coefs_1(end))]);
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


function YDATA = test_FMM(coefs_0, times, argsIn)
filename = {'./CamA/piv_ts.dfi', './CamB/piv_ts.dfi'};

YDATA = times*NaN;

for ii = 1:length(argsIn.StartEndInds)-1
    t_start = times(argsIn.StartEndInds(ii)+1);
    t_end = times(argsIn.StartEndInds(ii+1));

    if size(filename, 2)==1 %% Just be nice and load in one set of data
        im = dfireadvel(filename{1});
        u = im.cdata(:, :, 1);
        u(u == 0) = NaN; % Remove anomalous numbers
        u = flip(u, 1)';
        x = im.x;
        x = x(1, :)';
        timedata = flip(im.y(:, 1));

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
        timedata = new_t';
        x = new_x';
    else
        error('Length of filename wrong')
    end

    time_index = nearest_index(timedata, t_start):nearest_index(timedata, t_end);
    %timedata = timedata(time_index);
    u = u(:, time_index);

    Flow.U_flow = u;
    Flow.timestep = times(2)-times(1);
    Flow.x = x;

    Particle.StartLoc = argsIn.FunctionArgs.StartLoc(ii);
    Particle.C_d = coefs_0;
    Particle.rho_f = 910;
    Particle.Shape = 'circle';
    Particle.r = .150;
    [particle, ~] = FloatMotionModel(Flow, Particle, 'advanced');
    YDATA((argsIn.StartEndInds(ii)+1):argsIn.StartEndInds(ii+1)) = particle.x';

end
end
