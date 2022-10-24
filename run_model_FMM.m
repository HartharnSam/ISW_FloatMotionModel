%%run_model_FMM

clc; clearvars; close all; spinsstartup;
filename = {'./CamA/piv_ts.dfi', './CamB/piv_ts.dfi'};
t_start = 2; 

% Read in SPINS stuff
params = spins_params;
rho_1 = params.rho_0;
times = 0:length(dir('u.*'))-1;
time_index = nearest_index(times, t_start):max(times);
z_coord = -.01;
[~, x, u] = spins_hovmoller_z('u',z_coord, time_index, 'Plot', false);
times = times(time_index);

%% Parse & Run Model
Flow.U_flow = u;
Flow.timestep = times(2)-times(1);
Flow.x = x;
Flow.rho_0 = rho_1;

Particle.StartLoc = 4.77;
Particle.C_d = 17.77;
Particle.rho_f = 910;

[particle, fluid_u] = FloatMotionModel(Flow, Particle, 'advanced');
