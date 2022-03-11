%fit_ISWPTM
clc; clearvars; close all;
m_path = mpath;
%cd(m_path);
addpath(m_path);

%% Get the lab data
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
coefs_0 = [10];%, 910, locations(1)]; % [C_d, rho_f, x_start];
coefs_1 = lsqcurvefit(@test_APTM, coefs_0, times, locations);

[locations_fitted] = test_APTM(coefs_1, times);
fprintf('\n $C_d = $ %4.2f \n', coefs_1(1));

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


function x = test_APTM(coefs_0, times)
TypeParameter.Type  = 'lab';
TypeParameter.lab_fname = './CamA/piv_ts.dfi';
TypeParameter.Model = 'advanced';
TypeParameter.PlotType = 'line';
%[~, x] = advanced_PTM(TypeParameter, coefs_1(1), coefs_1(2), coefs_1(3));
[~, x] = advanced_PTM(TypeParameter, coefs_0(1), 910, 4.7921, false, 0, times(1));
x = x(round((times - times(1) + 1)*30))';

end
