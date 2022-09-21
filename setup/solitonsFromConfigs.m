function Psi = solitonsFromConfigs(simConfig)
	N = simConfig.N;
	Psi = cell(3, 1);
	for j = 1:3
		Psi{j} = zeros(N, N, N);
	end
	ctrs = simConfig.ctrs;
	epsilons = simConfig.epsilons;
	nSols = size(ctrs, 1);
	for i = 1:nSols
		if (simConfig.useExactProfiles)
			Sol = solitonExact(ctrs(i, :), simConfig.solIdxs(i), epsilons(i, :), simConfig);
		else
			Sol = solitonNodelessSi(ctrs(i, :), simConfig.r95s(i), epsilons(i, :), simConfig);
		end
		for j = 1:3
			Psi{j} = Psi{j} + Sol{j};
		end
	end
end