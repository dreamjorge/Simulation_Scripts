figure(1)
clf

% create panel
p = panel();
p.pack('h', {1/3 []})
p(1).pack(6, 2);

q = p(2);
for m = 1:6
	for n = 1:2
		
		% select axis - these two lines do the same thing (see
		% above)
% 		p(2, m, n).select();
		q(m, n).select();

		% prepare sample data
		data = randn(100, 1) * 0.4;
		
		% do stats
		stats = [];
		stats.source = source(m);
		stats.binrange = [-1 1];
		stats.xtick = [-0.8:0.4:0.8];
		stats.ytick = [0 20];
		stats.bincens = -0.9:0.2:0.9;
		stats.values = data;
		stats.freq = hist(data, stats.bincens);
		stats.percfreq = stats.freq / length(data) * 100;
		stats.percpeak = 30;
		
		% plot
		demopanel_minihist(stats, m == 6, n == 1);
		
  end
  
end