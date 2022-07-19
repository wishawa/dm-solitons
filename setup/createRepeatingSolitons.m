function Psi = createRepeatingSolitons(simConfig, createSolitonInGrid, repeatTimes)
	Psi = cell(1, 3);
	N = simConfig.N;
	dx = simConfig.dx;

	slin = (-N/2:N/2-1) * dx;
	[space1, space2, space3] = meshgrid(slin, slin, slin);
	Spaces = {space1, space2, space3};

	Lbox = simConfig.Lbox;
	gl = -repeatTimes:repeatTimes;
	[gx, gy, gz] = meshgrid(gl, gl, gl);
	go = {gx, gy, gz};
	for j = 1:3
		Psi{j} = zeros(N, N, N);
	end
	for i = 1:(length(gl) ^ 3)
		newSpaces = cell(1, 3);
		for j = 1:3
			newSpaces{j} = Spaces{j} + go{j}(i) * Lbox;
		end
		Soliton = createSolitonInGrid(newSpaces);
		for j = 1:3
			Psi{j} = Psi{j} + Soliton{j};
		end
	end
end