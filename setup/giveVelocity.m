function Out = giveVelocity(Soliton, Velocity, simConfig)
%myFun - Description
%
% Syntax: Out = giveVelocity(soliton)
%
% Long description
	N = simConfig.N;
	dx = simConfig.dx;

	slin = (-N/2:N/2-1) * dx;
	[space1, space2, space3] = meshgrid(slin, slin, slin);
	Spaces = {space1, space2, space3};

	Out = cell(1, 3);
	for j = 1:3
		Out{j} = Soliton{j};
		for k = 1:3
			Out{j} = Out{j} .* exp(1i * Spaces{k} .* Velocity(k));
		end
	end
end