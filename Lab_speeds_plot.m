clearvars; close all; clc;

filename = {'CamA/piv_ts.dfi', 'CamB/piv_ts.dfi'};

[u, x, timedata] = collate_piv_ts(filename);
load('CamA/ptv_tracks_compiled.mat')
c = 0.102;

particleR = .350/2;
timedata = timedata(1:max(length(ptv.n_timesteps), length(timedata)));

for ii = 1:length(timedata)
    uA(ii) = u(nearest_index(x, ptv.data{1}(ii)-particleR), ii);
    uB(ii) = u(nearest_index(x, ptv.data{1}(ii)+particleR), ii);

end
uA_sm = smooth(uA, 0.05);
uB_sm = smooth(uB, 0.05);
float_u = ptv.data{1}(1:max(length(ptv.n_timesteps), length(timedata)), 3);
tiledlayout(2, 1)
nexttile;
plot(timedata, uA/c, '-', 'Color',[1 0 0 .2]);
hold on
plot(timedata, uA_sm/c, '-', 'Color',[1 0 0]);

plot(timedata, uB/c, '-', 'Color',[0 0 1 .2]);
plot(timedata, uB_sm/c, '-', 'Color', [0 0 1]);

plot(timedata, float_u/c, 'k-')
ylim([-0.5 .5]);
yline(0, '-','Color', [1 1 1]*.3);
ylabel('$u/c_{isw}$', 'interpreter', 'latex')
legend('', 'Point A Fluid', '', 'Point B Fluid', 'Float',  'Location', 'best');
xticklabels([])
%title(['$L_f/\lambda = $ ', num2str(LfLambda)], 'interpreter', 'latex')

% Plot differences
nexttile;
plot(timedata, (float_u-uA_sm)/c, '-', 'Color',[1 0 0]);
hold on
plot(timedata, (float_u-uB_sm)/c, '-', 'Color', [0 0 1]);
ylim([-0.5 .5]);
yline(0, '-','Color', [1 1 1]*.3);
ylabel('$u_f - u(x) / c_{isw}$', 'interpreter', 'latex')
xlabel('t (s)')
xline([13 26], 'k')
xline([18 25], 'b')
figure_print_format(gcf);

exportgraphics(gcf, 'FloatFluidSpeeds.png')


function [u, x, timedata] = collate_piv_ts(filename)

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
end
