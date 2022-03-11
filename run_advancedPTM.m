clc; clearvars; close all;
%mpath = fileparts(mfilename('fullpath'));
m_path = mpath;
addpath([m_path, '../02_LabPlotting'])
cd([m_path, '../../02_Raw_data\CameraData\170222'])

TypeParameter.Type  = 'lab';
TypeParameter.lab_fname = './CamA/piv_ts.dfi';
TypeParameter.Model = 'advanced';
C_d = 10;
rho_f = 910;
x_start_loc = 4.77;
isPlot = false;

%% Load PTV data
load('./CamC/ptv_tracks.mat');
temp = ptv.data{1}(2, 3);
im = dfireadvel('./CamC/output_0000.dfi');
Grid = dfi_grid_read(im);
tempdata = interp1([1 Grid.nx], (Grid.x), temp);

u_0 = interp1([1 Grid.nx], (Grid.x), temp, 'linear', 'extrap') - Grid.x(1);

%%
[particle_t, particle_x, particle_u, u, x] = advanced_PTM(TypeParameter, C_d, rho_f, x_start_loc, isPlot, u_0);
pcolor(particle_t, x, u); shading flat; caxis([-.1 .1]);
newbluewhitered;
% colorbar;
hold on
plot(particle_t, particle_x);

plot_ptv_tracks('x', '030222');
