clearvars; close all; clc;

filename = {'CamA/piv_ts.dfi', 'CamB/piv_ts.dfi', 'CamC/piv_ts.dfi'};

[u, x, timedata] = collate_piv_ts(filename);
load('CamC/ptv_tracks_compiled.mat')
c = 0.107;

particleR = .350/2;
timedata = timedata(1:max(length(ptv.n_timesteps), length(timedata)));

uA = nan(1, length(timedata));
uB = uA;
for ii = 1:length(timedata)
    uA(ii) = u(nearest_index(x, ptv.data{1}(ii)-particleR), ii);
    uB(ii) = u(nearest_index(x, ptv.data{1}(ii)+particleR), ii);

end
uA_sm = smooth(uA, 0.05);
uB_sm = smooth(uB, 0.05);
float_u = ptv.data{1}(1:max(length(ptv.n_timesteps), length(timedata)), 3);
subaxis(2, 1, 1)

plot(timedata, uA/c, '-', 'Color',[1 0 0 .2]);
hold on
plot(timedata, uA_sm/c, '-', 'Color',[1 0 0]);

plot(timedata, uB/c, '-', 'Color',[0 0 1 .2]);
plot(timedata, uB_sm/c, '-', 'Color', [0 0 1]);

plot(timedata, float_u/c, 'k-')
ylim([-0.5 .5]);
yline(0, '-','Color', [1 1 1]*.3);
ylabel('$u/c_{isw}$', 'interpreter', 'latex')
legend('', ' Fluid', '', 'Point B Fluid', 'Float',  'Location', 'eastoutside');
xticklabels([])
xlim([0 70])
set(gca, 'Position', [0.1694 0.6098 0.5348 0.2757]);

%title(['$L_f/\lambda = $ ', num2str(LfLambda)], 'interpreter', 'latex')

% Plot differences
subaxis(2, 1, 2)
plot(timedata, (float_u-uA_sm)/c, '-', 'Color',[1 0 0]);
hold on
plot(timedata, (float_u-uB_sm)/c, '-', 'Color', [0 0 1]);
ylim([-0.5 .5]);
yline(0, '-','Color', [1 1 1]*.3);
ylabel('$u_f - u(x) / c_{isw}$', 'interpreter', 'latex')
xlabel('t (s)')
%xline([13 21 27 35 45 57 61], 'k')
xlim([0 70])
%xline([18 25], 'b')
set(gca, 'Position', [0.1694 0.1941 0.5348 0.2757]);

figure_print_format(gcf, 18)
fig = gcf;
fig.Units = 'centimeters';
fig.Position = [0 0 14 10.5];
%        exportgraphics(gcf,['../../04_Output/06_SurfaceFlow/BasicFlowFloatModel_', num2str(LfLambda), '.eps'], 'ContentType', 'vector')
%        exportgraphics(gcf,['../../04_Output/06_SurfaceFlow/BasicFlowFloatModel_', num2str(LfLambda), '.png'])
%exportgraphics(gcf,['../../04_Output/06_SurfaceFlow/FloatModels/BasicFlowFloatModel_', num2str(LfLambda), '.eps'], 'ContentType', 'vector')
dark_figure(gcf, [23 23 23])
export_fig(gcf,['FloatFluidSpeeds', '.png'], '-dpng')


function [u, x, timedata] = collate_piv_ts(filename)
cutoff = 8; % Amount to cut off the edges, which appears to be an effect on PIV images
for ii = 1:length(filename)
    im{ii} = dfireadvel(filename{ii});
    grid{ii} = dfi_grid_read(im{ii});
    im{ii}.cdata(:, [1:cutoff end-cutoff:end], :) = NaN;
    if ii == 1
        xmin = min(grid{ii}.x);
        xmax = max(grid{ii}.x);
        tmin = min(grid{ii}.y);
        tmax = max(grid{ii}.y);
    else
        xmin = min(xmin, min(grid{ii}.x));
        xmax = max(xmax, max(grid{ii}.x));
        tmin = min(tmin, min(grid{ii}.y));
        tmax = max(tmax, max(grid{ii}.y));
    end
end
% im1 = dfireadvel(filename{1});
% im2 = dfireadvel(filename{2});
% grid_1 = dfi_grid_read(im1);
% grid_2 = dfi_grid_read(im2);

% xmin= min(min(grid_1.x), min(grid_2.x));
% xmax= max(max(grid_1.x), max(grid_2.x));
% tmin= min(min(grid_1.y), min(grid_2.y));
% tmax= max(max(grid_1.y), max(grid_2.y));

new_x = xmax:grid{1}.dx:xmin;
new_t = tmin:grid{1}.dy:tmax;
[newX, newT] = meshgrid(new_x, new_t);
u = newX*NaN;
for ii = 1:length(filename)
    new_data = interp2(grid{ii}.X, grid{ii}.Y, im{ii}.cdata(:, :, 1), newX, newT);
    u(~isnan(new_data)) = new_data(~isnan(new_data));
end
u = u';
%new_1_data = interp2(grid_1.X, grid_1.Y, im1.cdata(:, :, 1), newX, newT);
%new_2_data = interp2(grid_2.X, grid_2.Y, im2.cdata(:, :, 1), newX, newT);
%u = new_1_data;
%u(~isnan(new_2_data)) = new_2_data(~isnan(new_2_data));
%u = u';
timedata = new_t';
x = new_x';
end
