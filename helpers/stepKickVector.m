function Psi = stepKickVector(Psi, dt, simConfig)
	%myFun - Description
	%
	% Syntax: Psi = halfKick(Psi, Vgrav, VsiScalar, dt, simConfig)
	%
	% Long description

	% [newPsi1, newPsi2, newPsi3] = arrayfun(@(psi1, psi2, psi3) updatePointPsi([psi1; psi2; psi3], dt, simConfig), Psi{1}, Psi{2}, Psi{3});
	% Psi = {newPsi1, newPsi2, newPsi3};

	OrigPsi = Psi;
	OrigPsiConj = conjVectorArray(OrigPsi);
	Psi = kickMainTerm(Psi, OrigPsi, dt / 2, simConfig);
	if (simConfig.doVectorCorrection)
		Psi = kickCorrectionTerm(Psi, OrigPsi, 1, dt, simConfig);
		Psi = kickCorrectionTerm(Psi, OrigPsiConj, -1, dt, simConfig);
	end
	Psi = kickMainTerm(Psi, OrigPsi, dt / 2, simConfig);
	% if (simConfig.doVectorCorrection)
	% 	Psi = kickCorrectionTerm(Psi, OrigPsiConj, -1, dt / 2, simConfig);
	% 	Psi = kickCorrectionTerm(Psi, OrigPsi, 1, dt / 2, simConfig);
	% end
end
function NewPsi = kickCorrectionTerm(Psi, PsiCc, sgn, dt, simConfig)
	PsiCcSq = dotVectorArray(PsiCc, PsiCc);
	PsiCcSqAbs = abs(PsiCcSq);
	PsiCcSqAbsSq = PsiCcSqAbs .^ 2;
	NewPsi = Psi;
	MCoef = (exp(sgn * dt^2 * simConfig.siCoef^2 * 0.5 * PsiCcSqAbsSq) - 1) ./ PsiCcSq;
	MCoef(PsiCcSqAbs < 1E-12) = 1;
	% MCoef(isnan(MCoef)) = 0;
	MCoef = MCoef .* dotVectorArray(PsiCc, Psi);
	for j = 1:3
		NewPsi{j} = NewPsi{j} + MCoef .* PsiCc{j};
	end
	% for j = 1:3
	% 	for k = 1:3
	% 		NewPsi{j} = NewPsi{j} + (MCoef .* PsiCc{j} .* PsiCc{k}) .* Psi{k};
	% 	end
	% end
end
function NewPsi = kickMainTerm(Psi, PsiForOp, dt, simConfig)
	rhoMul = -1i * dt * simConfig.siCoef;
	Rho = getRho(PsiForOp, simConfig);
	MCoef = (exp(Rho * rhoMul) - 1) ./ Rho;
	MCoef(Rho < 1E-12) = 1;
	NewPsi = Psi;
	MCoef = MCoef .* dotVectorArray(PsiForOp, Psi);
	for j = 1:3
		NewPsi{j} = NewPsi{j} + MCoef .* conj(PsiForOp{j});
	end
	% for j = 1:3
	% 	for k = 1:3
	% 		NewPsi{j} = NewPsi{j} + (MCoef .* conj(PsiForOp{j}) .* PsiForOp{k}) .* Psi{k};
	% 	end
	% end
end