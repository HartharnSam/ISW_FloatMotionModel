clc; clearvars; close all; 
addpath('./')

addpath(genpath('../ISWLab_Toolkit'));
cd('../../02_Raw_data/CameraData/170222');
h_1 = 0.07;
h_2 = 0.23;
h_ice = 0.03;
h_pyc = 0.015;
wave_amp = 0.077;
rho_1 = 1028;
rho_2 = 1049;

[wavelength, DJL] = calc_DJL(h_1, h_2, h_pyc, wave_amp, rho_1, rho_2);
%% Plot out DJL results and identify ice lower level
z_ind = nearest_index(DJL.z, -h_ice);
clf; 
pcolor(DJL.x, DJL.z, DJL.density); shading flat; caxis([.99 1.01]); cmocean('dense')
hold on
contour(DJL.x, DJL.z, DJL.u, 'k-'); 
yline(DJL.z(z_ind), '--r');

%% 
u_max = max(DJL.u(z_ind, :));
wave_center = 0;
c_w = DJL.WaveC;
negative_amplitude = u_max*.4;
negative_wave_center = wave_center - wavelength/2;
c_neg_w = c_w;
wavelength_negative = wavelength *.1;
zero_offset = 1;

%0.534011663726336,0.00967525983580827,-1.67876599222030,0.145593645787183,0.999999999999970]

input_parameters = [u_max, wave_center, c_w, wavelength/4, negative_amplitude,...
    negative_wave_center, c_neg_w, wavelength_negative, zero_offset];
Times = [0:60];
[XDATA(:, :, 1) XDATA(:, :, 2)] = meshgrid(DJL.x, Times);

FittedData = ISWFitFunction(input_parameters, XDATA);
figure;
pcolor(DJL.x, Times, FittedData); caxis([-.065 .065])
shading flat;
cmocean('balance');
