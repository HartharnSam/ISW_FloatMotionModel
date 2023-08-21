%%run_lab_FMM - Parsing Script for lab Data into the FMM
%u = FitLabFlow(TypeParameter.lab_fname, 1);
% makes some assumptions about the input piv_ts:
%   Start at same time, and have the same timesteps
%   That it's okay for the second input camera to overwrite (rather than a
%   complicated merge with) CamA for region of overlap
%
clc; clearvars; close all; digiflowstartup;
filename = {'./CamA/piv_ts.dfi', './CamB/piv_ts.dfi'};
isFillMissing = false; 
t_start = 2; 

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

Particle.StartLoc = 4.77;
Particle.C_d = 17.77;
Particle.rho_f = 910;

[particle, fluid_u] = FloatMotionModel(Flow, Particle, 'advanced');
