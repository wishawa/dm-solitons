function totalMass = getTotalMass(Rho, simConfig)
%myFun - Description
%
% Syntax: totalMass = getTotalMass(Rho)
%
% Long description
totalMass = sum(Rho, 'all') * simConfig.dx^3;
end