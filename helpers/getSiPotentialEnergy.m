function EVsi = getSiPotentialEnergy(Psi, simConfig)
%myFun - Description
%
% Syntax: EVsi = getSiPotentialEnergy(Psi, simConfig)
%
% Long description
	
N = simConfig.N;
PsiSq = zeros(N, N, N);
Rho = zeros(N, N, N);
for j = 1:3
	PsiSq = PsiSq + Psi{j}.^2;
	Rho = Rho + abs(Psi{j}).^2;
end
EVsi = 0.25 * simConfig.lambda * sum((2 * (Rho.^2) + abs(PsiSq).^2), 'all') * simConfig.dx^3;
end