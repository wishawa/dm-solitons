function ET = getKineticEnergy(Psi, kGrids, simConfig)
%myFun - Description
%
% Syntax: T = getKineticEnergy(Psi)
%
% Long description

T = 0;
for j = 1:3
	FourierPsi = fftn(Psi{j});
	for k = 1:3
		partialJK = ifftn(1i * kGrids{k} .* FourierPsi);
		T = T + sum(abs(partialJK).^2, 'all');
	end
end
ET = 0.5 * T * simConfig.dx^3;
end