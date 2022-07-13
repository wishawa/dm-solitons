function ET = getKineticEnergy(Psi, kGrids, simConsts)
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
ET = T * simConsts.dx^3 / (2 * simConsts.m_per_hbar^2);
end