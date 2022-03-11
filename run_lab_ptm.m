%%run_lab_PTM;
%u = FitLabFlow(TypeParameter.lab_fname, 1);
filename = 'piv_ts.dfi';

im = dfireadvel(TypeParameter.lab_fname);
u = im.cdata(:, :, 1);
u(u == 0) = NaN; % Remove anomalous numbers
%u = fillmissing(u, 'linear', 'EndValues', 'none')'; %Re-fill those values
u = flip(u, 1)';

xi = im.x;
xi = xi(1, :)';
times = flip(im.y(:, 1));
time_index = nearest_index(times, t_start):length(times);
times = times(time_index);
u = u(:, time_index);