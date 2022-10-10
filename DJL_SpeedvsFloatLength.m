% Plots calculated float motion according to DJL waves and the PTM model

LfLambda = [0.01:.01:1]; %Float length/Wavelength
cf_cw = [0.400608652000000	0.399165065000000	0.396830606000000	0.394014650000000	0.390025406000000	0.385186724000000	0.379528388000000	0.373842511000000	0.366756181000000	0.359002609000000	0.351614298000000	0.342826533000000	0.333631409000000	0.324136553000000	0.315530680000000	0.305760386000000	0.295991324000000	0.287380529000000	0.277845819000000	0.268533737000000	0.259494704000000	0.251719354000000	0.243287554000000	0.235204738000000	0.228319120000000	0.220909975000000	0.213853440000000	0.207140869000000	0.201453621000000	0.195357676000000	0.189567572000000	0.184665156000000	0.179410692000000	0.174417388000000	0.169670150000000	0.165645130000000	0.161323311000000	0.157207151000000	0.153710782000000	0.149949273000000	0.146359044000000	0.142929620000000	0.140008327000000	0.136856525000000	0.133839094000000	0.131263290000000	0.128478415000000	0.125806412000000	0.123240839000000	0.121044809000000	0.118664218000000	0.116373819000000	0.114409664000000	0.112276560000000	0.110220410000000	0.108237218000000	0.106532600000000	0.104677212000000	0.102884636000000	0.101341468000000	0.0996592730000000	0.0980314960000000	0.0964555730000000	0.0950963340000000	0.0936119050000000	0.0921727580000000	0.0907768760000000	0.0895708710000000	0.0882516050000000	0.0869703990000000	0.0858621930000000	0.0846485510000000	0.0834685550000000	0.0823208340000000	0.0813266800000000	0.0802364380000000	0.0791749060000000	0.0782545220000000	0.0772442280000000	0.0762595820000000	0.0752996230000000	0.0744663130000000	0.0735505330000000	0.0726569230000000	0.0718805740000000	0.0710267120000000	0.0701928320000000	0.0693782440000000	0.0686698270000000	0.0678899040000000	0.0671274480000000	0.0664639050000000	0.0657328850000000	0.0650177280000000	0.0643179270000000	0.0637083760000000	0.0630362640000000	0.0623781520000000	0.0618045670000000	0.0611717370000000];
% max float speed / c_isw

figure
plot(LfLambda, cf_cw);
xline(0.04, 'r:')
xlabel('$L_f / \lambda$', 'interpreter', 'latex')
ylabel('$Max(c_f/c_{isw})$', 'interpreter', 'latex');

hold on
% Add in the data from Carr et al 2019 experiments (pulled from plot
% digitiser, figure 3b)
x = [0.7616055846422339, 0.7895287958115184, 1.0226876090750436, 0.5465968586387434,0.5396160558464224];
z = [0.35395348837209306, 0.3074418604651163, 0.15209302325581392,  0.4451162790697675,  0.4534883720930233];

plot(x, z, 'xr')