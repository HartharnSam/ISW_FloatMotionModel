% Plots calculated float motion according to DJL waves and the PTM model
clc; clearvars; close all;

%% Load in data
load max_uf_c
labWavelengths = [2.052, 2.523, 1.572, 2.072, 1.546, 1.469, 1.456, 2.299, 1.705, 1.955, 1.561, 2.863, 1.835];
labData = [0.064475, 0.066037, 0.049031, 0.071453, 0.061702, 0.033281, 0.022739, 0.045763, NaN, 0.075193, 0.045461, 0.0526, 0.0166];
labC = [0.128, 0.125, 0.118, 0.102, 0.109, 0.111, 0.12, 0.114, 0.11, 0.106, 0.107 , 0.136, 0.128];
labData = labData./labC;
labLfLambda = [0.0604607, 0.05358, 0.06511727, 0.1523, 0.1999, 0.2163, 0.2256, 0.17105, 0.1844, 0.1393, 0.2119, 0.4192, 0.654];

% Sort colormaps
figure
cmap = plasma;
%max_uf_c.lambda = max_uf_c.amp; 

wavelengths = [min(min(max_uf_c.lambda), min(labWavelengths)), max(max(max_uf_c.lambda), max(labWavelengths))];
wavelengths = linspace(wavelengths(1), wavelengths(2), length(cmap));

for ii = 1:length(max_uf_c.lambda)
    cmap_ind = nearest_index(wavelengths, max_uf_c.lambda(ii));
    semilogx(max_uf_c.LfLambda, max_uf_c.data(:, ii), '-', 'Color', cmap(cmap_ind, :), 'DisplayName', strcat('$\lambda = ',num2str(max_uf_c.lambda(ii)'), 'm$'));
    hold on

end
xlabel('$L_f / \lambda$', 'interpreter', 'latex')
ylabel('$Max(u_f/c_{isw})$', 'interpreter', 'latex');
%legend('interpreter', 'latex', 'AutoUpdate','off');

hold on
for ii = 1:length(labData)
    cmap_ind = nearest_index(wavelengths, labWavelengths(ii));
    plot(labLfLambda(ii), labData(ii), 'x', 'Color', cmap(cmap_ind, :))

end

% Add in the data from Carr et al 2019 experiments (pulled from plot
% digitiser, figure 3b)
%x = [0.7616055846422339, 0.7895287958115184, 1.0226876090750436, 0.5465968586387434,0.5396160558464224];
%z = [0.35395348837209306, 0.3074418604651163, 0.15209302325581392,  0.4451162790697675,  0.4534883720930233];

%plot(x, z, 'xr');

figure_print_format(gcf, 18)
fig = gcf;
fig.Units = 'centimeters';
fig.Position = [0 0 14 10.5];
axis tight
dark_figure(gcf, [23 23 23])
export_fig(gcf,['LfvsUmax.png'], '-dpng')
