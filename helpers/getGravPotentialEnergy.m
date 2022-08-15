function EVgrav = getGravPotentialEnergy(Vgrav, Rho)
%myFun - Description
%
% Syntax: Vgrav = getGravPotentialEnergy(Vgrav, Rho, simConfig)
%
% Long description
EVgrav = 0.5 * sum(Vgrav .* Rho, 'all') * simConfig.dx^3;
end