function Psi = stepKickVector(Psi, dt, simConsts)
	%myFun - Description
	%
	% Syntax: Psi = halfKick(Psi, Vgrav, VsiScalar, dt, simConsts)
	%
	% Long description

	% [newPsi1, newPsi2, newPsi3] = arrayfun(@(psi1, psi2, psi3) updatePointPsi([psi1; psi2; psi3], dt, simConsts), Psi{1}, Psi{2}, Psi{3});
	% Psi = {newPsi1, newPsi2, newPsi3};

	Psi = kickCorrectionTerm(Psi, Psi, 1, dt, simConsts);
	Psi = kickCorrectionTerm(Psi, conjVectorArray(Psi), -1, dt, simConsts);
	Psi = kickMainTerm(Psi, dt, simConsts);
	% Psi = kickCorrectionTerm(Psi, conjVectorArray(Psi), -1, dt / 2, simConsts);
	% Psi = kickCorrectionTerm(Psi, Psi, 1, dt / 2, simConsts);

end
function NewPsi = kickCorrectionTerm(Psi, PsiCc, sgn, dt, simConsts)
	PsiCcSq = dotVectorArray(PsiCc);
	PsiCcSqDag = conj(PsiCcSq);
	NewPsi = Psi;
	MCoef = (exp(sgn * dt^2 * simConsts.siCoef^2 * 0.5 * PsiCcSqDag .* PsiCcSq) - 1) ./ PsiCcSq;
	for j = 1:3
		for k = 1:3
			NewPsi{j} = NewPsi{j} + (MCoef .* PsiCc{j} .* PsiCc{k}) .* Psi{k};
		end
	end
end
function NewPsi = kickMainTerm(Psi, dt, simConsts)
	rhoMul = -1i * dt * simConsts.siCoef;
	Rho = getRho(Psi, simConsts);
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
	psiUpdated = expm(simConsts.siCoef^2 * kickDt^2 * 0.5 * dHdt) * psi;
	o1 = psiUpdated(1);
	o2 = psiUpdated(2);
	o3 = psiUpdated(3);
end