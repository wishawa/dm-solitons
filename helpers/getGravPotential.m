function Vgrav = getGravPotential(Rho, rhobar, kSqNonzero, simConfig)
%myFun - Description
%
% Syntax: Vscalar = getScalarPotential(Rho, rhobar, simConfig)
%
% Long description

% Gravitational Potential
Vgrav = -4 * pi * simConfig.G * simConfig.m_per_hbar * (Rho - rhobar);
Vgrav = fftn(Vgrav);
Vgrav = Vgrav ./ kSqNonzero;
Vgrav = ifftn(Vgrav);

end