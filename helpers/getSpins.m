function Spins = getSpins(Psi, simConfig, normalized)
%spinRaw - Compute a 3x1 cells array of double arrays containing the spin in each component.
%
% Syntax: Spins = spinRaw(Psi)
%
% Long description
	
Spins = cell(1, 3);
for j = 1:3
	PsiB = Psi{rem(j, 3) + 1};
	PsiC = Psi{rem(j + 1, 3) + 1};
	if normalized
		normConst = 1.;
	else
		normConst = simConfig.dx^3;
	end
	Spins{j} = imag(conj(PsiB) .* PsiC - PsiB .* conj(PsiC)) * normConst;
end
end