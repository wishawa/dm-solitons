function EVgrav = getGravPotentialEnergy(Vgrav, Rho, simConfig)
%myFun - Description
%
% Syntax: Vgrav = getGravPotentialEnergy(Vgrav, Rho, simConfig)
%
% Long description
EVgrav = sum(Vgrav .* Rho, 'all') * simConfig.dx^3 / (2 * simConfig.m_per_hbar);
end