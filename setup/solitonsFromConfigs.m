function Psi = solitonsFromConfigs(simConfig)
	N = simConfig.N;
	Psi = cell(3, 1);
	for j = 1:3
		Psi{j} = zeros(N, N, N);
	end
	ctrs = simConfig.ctrs;
	r95s = simConfig.r95s;
	% solIdxs = simConfig.solIdxs;
	epsilons = simConfig.epsilons;
	nSols = size(ctrs, 1);
	for i = 1:nSols
		Sol = solitonNodelessSi(ctrs(i, :), r95s(i), epsilons(i, :), simConfig);
		% Sol = solitonExact(ctrs(i, :), solIdxs(i), epsilons(i, :), simConfig);
		for j = 1:3
			Psi{j} = Psi{j} + Sol{j};
		end
	end
end