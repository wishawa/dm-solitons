function Psi = createRepeatingSolitons(simConsts, createSolitonsInGrid)
	Psi = cell(1, 3);
	N = simConsts.N;
	dx = simConsts.dx;

	slin = (-N/2:N/2-1) * dx;
	[space1, space2, space3] = meshgrid(slin, slin, slin);
	Spaces = {space1, space2, space3};

	Lbox = simConsts.Lbox;
	gl = -1:1;
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
		Solitons = createSolitonsInGrid(newSpaces, simConsts);
		for k = 1:length(Solitons)
			sk = Solitons{k};
			for j = 1:3
				Psi{j} = Psi{j} + sk{j};
			end
		end
	end
end