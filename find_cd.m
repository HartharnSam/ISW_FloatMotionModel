%FIT_CD
% Fit C_d 
clc; clearvars; close all; digiflowstartup;

m_path = mpath;
%cd(m_path);
addpath(m_path);

%% Get the lab Particle Data
load('./CamC/ptv_tracks.mat', 'ptv'); % Load in the data
im = dfireadvel('CamC/output_0000.dfi'); % Load in a frame to get WCS adjustment
Grid = dfi_grid_read(im);
ii = 1;
times = (0:ptv.n_timesteps-1)+1; % Time indices

%TODO Collate the tracks as if one
for i = 1:ptv.n_particles
    if ~isempty(ptv.data{i})
        locations_tmp = interp1([1 Grid.nx], (Grid.x), ptv.data{i}(:, 32));
        locations_tmp = locations_tmp(~isnan(locations_tmp))'; % Remove nans
        locations{ii} = locations_tmp;
        times_all{ii} = times(~isnan(locations_tmp))/30;

        %% Make first guess at least squares fit
        coefs_0 = 10;%, 910, locations(1)]; % [C_d, rho_f, x_start];
        t_starts(ii) = nearest_index(times_all{ii}, t_start);
        t_ends(ii) = nearest_index(times_all{ii}, t_end);
        times_all{ii} = times_all{ii}(t_starts(ii):t_ends(ii));
        locations{ii} = locations{ii}(t_starts(ii):t_ends(ii));
        loc_starts(ii) = locations{ii}(1);
        ii = ii+1;

    end
end
function_args.tstarts = t_starts;
function_args.t_ends = t_ends;
function_args.StartLoc = loc_starts;
[coefs_1,~,~,~,~,~,~, argsIn] = lsqcurvefit_ul(@test_FMM, coefs_0, times_all, locations, function_args); % TODO: Make it not the final point
%%
num_times = 0;
for ii = 1:length(times_all)
    num_times = num_times+length(times_all{ii});
end
xdata = NaN(num_times, 1); ydata = xdata;
end_index = 0;
for ii = 1:length(times_all)
    indexes = end_index(ii)+(1:length(times_all{ii}));
    end_index(ii+1) = indexes(end);
    xdata(indexes) = times_all{ii};
    ydata(indexes) = locations{ii};
end
times_2 = xdata; 
locations_2 = ydata;
[locations_fitted] = test_FMM(coefs_1, times_2, argsIn);
fprintf('\n C_d =  %4.2f \n', coefs_1(1));

%% Plot the data
figure

plot(locations_2, times_2, 'kx');  % Plot the lab data
hold on
plot(locations_fitted, times_2, 'b-')

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
    timedata = timedata(time_index);
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
