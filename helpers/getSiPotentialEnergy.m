function EVsi = getSiPotentialEnergy(Psi, simConsts)
%myFun - Description
%
% Syntax: EVsi = getSiPotentialEnergy(Psi, simConsts)
%
% Long description
	
coef = simConsts.siCoef / (2 * simConsts.m_per_hbar);
N = simConsts.N;
PsiSq = zeros(N, N, N);
Rho = zeros(N, N, N);
for j = 1:3
	PsiSq = PsiSq + Psi{j}.^2;
	Rho = Rho + abs(Psi{j}).^2;
end
EVsi = coef * sum((2 * (Rho.^2) + abs(PsiSq).^2), 'all');
end