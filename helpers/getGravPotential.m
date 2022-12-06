function Vgrav = getGravPotential(Rho, kSqNonzero)
%myFun - Description
%
% Syntax: Vscalar = getScalarPotential(Rho, rhobar, simConfig)
%
% Long description

% Gravitational Potential
Vgrav = fftn(-0.5 * Rho);
Vgrav = Vgrav ./ kSqNonzero;
Vgrav = ifftn(Vgrav);

end