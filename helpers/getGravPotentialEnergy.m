function EVgrav = getGravPotentialEnergy(Vgrav, Rho, simConsts)
%myFun - Description
%
% Syntax: Vgrav = getGravPotentialEnergy(Vgrav, Rho, simConsts)
%
% Long description
EVgrav = sum(Vgrav .* Rho, 'all') * simConsts.dx^3 / (2 * simConsts.m_per_hbar);
end