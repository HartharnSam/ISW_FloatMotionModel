close all; clc; clearvars;
params = spins_params;
x_start_loc = 4.68;
[t, x] = advanced_PTM(.2, 68, x_start_loc, true);
figure
plot(t, x);

filename = 'C:\Users\samha\OneDrive - Newcastle University\02_PhD_Project\02_Ice_Covered_Waters\02_Raw_data\CameraData\150921\CamC\particles.dfd';
tracked_particles = [2];
p = dfdread(filename);
xexp = nan(length(1:1500), length(tracked_particles));
for i = 1:1500
    for j = 1:length(tracked_particles)
        particle = tracked_particles(j)
        try
            xexp(i, j) = p.Data{i}(particle, 1);
        end
    end
    
    texp(i) = i/30+25;
end
hold on
for j =  1:length(tracked_particles)
    plot(texp, smooth(xexp(:, j), 150, 'lowess'))
    hold on
    legend
end

figure(1)
hold on
pi = plot(smooth(xexp(:, j), 150, 'lowess'), texp, 'r-')


%print('../../../05_Output/figure.png', '-dpng')