function NewPsi = stepKickVector(Psi, Rho, kickDt, simConsts)
%myFun - Description
%
% Syntax: Psi = halfKick(Psi, Vgrav, VsiScalar, dt, simConsts)
%
% Long description

% dHdt = cell(3, 3);
% dHdtCoef = kickDt^2 * simConsts.siCoef^2 / 2;
% PsiDag = {conj(Psi{1}), conj(Psi{2}), conj(Psi{3})};
% PsiSq = squareVectorArrays(Psi);
% PsiDagSq = conj(PsiSq);
% for j = 1:3
% 	for k = 1:3
% 		dHdt{j, k} = dHdtCoef * (PsiDagSq .* Psi{j} .* Psi{k} - PsiDag{j} .* PsiDag{k} .* PsiSq);
% 	end
% end

[newPsi1, newPsi2, newPsi3] = arrayfun(@(psi1, psi2, psi3) updatePointPsi([psi1; psi2; psi3], kickDt, simConsts), Psi{1}, Psi{2}, Psi{3});
Psi = {newPsi1, newPsi2, newPsi3};

rhoMul = -1i * kickDt * simConsts.siCoef;
MCoef = (exp(Rho * rhoMul) - 1) ./ Rho;
NewPsi = Psi;
for j = 1:3
	for k = 1:3
		NewPsi{j} = NewPsi{j} + (MCoef .* conj(Psi{j}) .* Psi{k}) .* Psi{k};
	end
end

end
function [o1, o2, o3] = updatePointPsi(psi, kickDt, simConsts)
	psiSq = dot(psi, psi);
	psiSqDag = conj(psiSq);

	dHdt = psi * transpose(psi) * psiSqDag - conj(psi) * transpose(conj(psi)) * psiSq;
	psiUpdated = expm(simConsts.siCoef^2 * kickDt^2 * -0.5 * dHdt) * psi;
	o1 = psiUpdated(1);
	o2 = psiUpdated(2);
	o3 = psiUpdated(3);
end