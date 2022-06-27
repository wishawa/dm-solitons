function Vscalar = getSiScalarPotential(Rho, simConsts)
%myFun - Description
%
% Syntax: Vscalar = getScalarPotential(Rho, rhobar, simConsts)
%
% Long description

% Scalar SI Potential
Vscalar = simConsts.siCoef * 2 * Rho;

end