function Vgrav = getGravPotential(Rho, rhobar, kSqNonzero, simConsts)
%myFun - Description
%
% Syntax: Vscalar = getScalarPotential(Rho, rhobar, simConsts)
%
% Long description

% Gravitational Potential
Vgrav = -4 * pi * simConsts.G * simConsts.m_per_hbar * (Rho - rhobar);
Vgrav = fftn(Vgrav);
Vgrav = Vgrav ./ kSqNonzero;
Vgrav = ifftn(Vgrav);

end