function Psi = solitonsFromConfigs(simConfig)
	N = simConfig.N;
	Psi = cell(3, 1);
	for j = 1:3
		Psi{j} = zeros(N, N, N);
	end
	positions = simConfig.positions;
	sizes = simConfig.sizes;
	spins = simConfig.spins;
	nSols = size(simConfig.positions, 1);
	for i = 1:nSols
		Sol = solitonNodelessSi(positions(i, :), sizes(i), spins(i, :), simConfig);
		for j = 1:3
			Psi{j} = Psi{j} + Sol{j};
		end
	end
end