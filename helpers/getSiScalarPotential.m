function Vscalar = getSiScalarPotential(Rho, simConfig)
%myFun - Description
%
% Syntax: Vscalar = getScalarPotential(Rho, rhobar, simConfig)
%
% Long description

% Scalar SI Potential
Vscalar = simConfig.lambda * Rho;

end