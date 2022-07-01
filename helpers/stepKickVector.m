function NewPsi = stepKickVector(Psi, PsiForVsi, Rho, kickDt, simConsts)
%myFun - Description
%
% Syntax: Psi = halfKick(Psi, Vgrav, VsiScalar, dt, simConsts)
%
% Long description

NewPsi = Psi;

rhoMul = -1i * kickDt * simConsts.siCoef;
MCoef = (exp(Rho * rhoMul) - 1) ./ Rho;
for j = 1:3
	for k = 1:3
		NewPsi{j} = NewPsi{j} + (MCoef .* conj(PsiForVsi{j}) .* PsiForVsi{k}) .* Psi{k};
	end
end
end