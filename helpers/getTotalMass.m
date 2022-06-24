function totalMass = getTotalMass(Rho, simConsts)
%myFun - Description
%
% Syntax: totalMass = getTotalMass(Rho)
%
% Long description
totalMass = sum(Rho, 'all') * simConsts.dx^3;
end