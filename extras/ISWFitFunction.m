function F = ISWFitFunction(a,XDATA)
%ISWFITFUNCTION - Generalised function for an ISW velocity field
% Gaussian profile in x with the location parameter changing linearly in y
% (time)
%
% Inputs:
%    a - 9x1 matrix of coefficients [amplitude wave_centre(1) c_w wavelength...
%           negative_amplitude negative_wave_center c_neg_w wavelength_negative zero_offset];
%    XDATA - Matrix of x, z from meshgrid as [XDATA(:, :, 1), XDATA(:, :,
%    2)]
%    = meshgrid(x, z);
%
% Outputs:
%    F - Data calculated from X-Y plane in XDATA and the coefficients in a
%    output2 - Description
%
% Example:
%   [time,x, data ] = spins_hovmoller_z('u', -.03, [10:60]);
%   [X TIME] = meshgrid(x, time);
%   XDATA(:, :, 1) = X;
%   XDATA(:, :, 2) = TIME;
%   a =[0.0522, .2552, 1.01, .09, .1]; % Fit parameters for 

%   DATA = ISWFitFunction(a, XDATA);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: FitFlow
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 02-Feb-2022; Last revision: 02-Feb-2022
% MATLAB Version: 9.10.0.1602886 (R2021a)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
xdata = XDATA(:, :, 1);
ydata = XDATA(:, :, 2);

%F = a(1)*exp(-((xdata-(a(2) + a(3)*ydata)).^2/(2*a(4)^2))) - ...
%    a(5)*exp(-((xdata-(a(6) + a(7)*ydata)).^2/(2*a(8)^2))) - a(9);

F = a(1)*sech(-((xdata-(a(2) + a(3)*ydata))/(a(4)))).^2 - ...
    a(5)*sech(-((xdata-(a(6) + a(7)*ydata))/(a(8)))).^2 - a(9);


end
