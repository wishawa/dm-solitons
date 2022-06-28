function EVsi = getSiPotentialEnergy(Psi, simConsts)
%myFun - Description
%
% Syntax: EVsi = getSiPotentialEnergy(Psi, simConsts)
%
% Long description
	
N = simConsts.N;
PsiSq = zeros(N, N, N);
Rho = zeros(N, N, N);
for j = 1:3
	PsiSq = PsiSq + Psi{j}.^2;
	Rho = Rho + abs(Psi{j}).^2;
end
EVsi = simConsts.siCoef * sum((2 * (Rho.^2) + abs(PsiSq).^2), 'all') * simConsts.dx^3 / (2 * simConsts.m_per_hbar);
end