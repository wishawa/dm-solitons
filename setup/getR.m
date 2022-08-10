function R = getR(ctr, simConfig)
%myFun - Description
%
% Syntax: R = getR(Spaces, center)
%
% Long description
	N = simConfig.N;
	dx = simConfig.dx;
	Lbox = simConfig.Lbox;

	slin = (-N/2:N/2-1) * dx;
	[space1, space2, space3] = meshgrid(slin, slin, slin);
	Spaces = {space1, space2, space3};

	R = zeros(N, N, N);
	for j = 1:3
		di = abs(Spaces{j} - ctr(j));
		di = min(di, Lbox - di);
		R = R + di.^2;
	end
	R = sqrt(R);
end