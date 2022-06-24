function NewPsi = fullDrift(Psi, kSq, dt, simConsts)
%myFun - Description
%
% Syntax: Psi = fullDrift(Psi, kSq, dt, simConsts)
%
% Long description
driftCoef = -1i * dt / (2 * simConsts.m_per_hbar) * kSq;
NewPsi = cell(1, 3);
for j = 1:3
	FourierPsi = fftn(Psi{j});
	FourierPsi = exp(driftCoef) .* FourierPsi;
	NewPsi{j} = ifftn(FourierPsi);
end
	
end