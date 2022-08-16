function Spins = getSpins(Psi, simConfig)
%spinRaw - Compute a 3x1 cells array of double arrays containing the spin in each component.
%
% Syntax: Spins = spinRaw(Psi)
%
% Long description
	
Spins = cell(1, 3);
for j = 1:3
	PsiB = Psi{rem(j, 3) + 1};
	PsiC = Psi{rem(j + 1, 3) + 1};
	Spins{j} = imag(conj(PsiB) .* PsiC - PsiB .* conj(PsiC)) * simConfig.dx^3;
end
end