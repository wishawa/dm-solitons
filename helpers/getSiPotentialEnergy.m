function EVsi = getSiPotentialEnergy(Psi, simConsts)
%myFun - Description
%
% Syntax: EVsi = getSiPotentialEnergy(Psi, simConsts)
%
% Long description
	
EVsi = 0;
coef = simConsts.siCoef / (2 * simConsts.m_per_hbar);
for j = 1:3
	term1 = sum(abs(Psi{j}).^4, 'all') * 2;
	term2 = sum(abs(Psi{j}.^2).^2, 'all');
	EVsi = EVsi + (term1 + term2);
end
EVsi = EVsi * coef * simConsts.dx^3;
end