function NewPsi = stepKickScalar(Psi, VScalar, kickDt)
%myFun - Description
%
% Syntax: Psi = halfKick(Psi, Vgrav, VsiScalar, dt, simConfig)
%
% Long description

NewPsi = Psi;

for j = 1:3
	NewPsi{j} = exp(-1i * kickDt * VScalar) .* NewPsi{j};
end
end