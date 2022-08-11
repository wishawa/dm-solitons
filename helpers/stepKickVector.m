function Psi = stepKickVector(Psi, dt, simConfig)
	%myFun - Description
	%
	% Syntax: Psi = halfKick(Psi, Vgrav, VsiScalar, dt, simConfig)
	%
	% Long description

	% [newPsi1, newPsi2, newPsi3] = arrayfun(@(psi1, psi2, psi3) updatePointPsi(psi1, psi2, psi3, dt, simConfig), Psi{1}, Psi{2}, Psi{3});
	% Psi = {newPsi1, newPsi2, newPsi3};
	% N = simConfig.N;
	% Psi1 = Psi{1};
	% Psi2 = Psi{2};
	% Psi3 = Psi{3};
	% for l = 1:(N^3)
	% 	% [p1, p2, p3] = updatePointPsi(Psi1(l), Psi2(l), Psi3(l), dt, simConfig);
	% 	psi = [Psi1(l); Psi2(l); Psi3(l)];
	% 	psiSq = dot(psi, psi);
	% 	psiSqDag = conj(psiSq);
	% 	psiDag = conj(psi);
	% 	matr = (psi * transpose(psi) * psiSqDag - psiDag * transpose(psiDag) * psiSq) * (simConfig.siCoef^2 * dt^2 * -0.5) + (psiDag * transpose(psi)) * (simConfig.siCoef * dt * -1i);
	% 	psiUpdated = expm(matr) * psi;
	% 	Psi1(l) = psiUpdated(1);
	% 	Psi2(l) = psiUpdated(2);
	% 	Psi3(l) = psiUpdated(3);
	% end

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
function EigVec = evecFromLambda(A11, A12, A13, A21, A22, A23, A33, lambda)
	M11 = A11 - lambda;
	M22 = A22 - lambda;
	% M33 = A33 - lambda;
	Can1 = cross([M11; A12; A13], [A12; M22; A23]);
	EigVec = Can1;
	% Can2 = cross([M11; A12; A13], [A12; M22; A23]);
	% Can3 = cross([M11; A12; A13], [A12; M22; A23]);
end
function expLattice(Psi, simConfig)
	PsiConj = {conj(Psi{1}), conj(Psi{2}), conj(Psi{3})};
	siCoef = simConfig.siCoef;
	MCoef = -0.5 * siCoef^2 * dt^2;

	PsiSq = dotVectorArray(Psi, Psi);
	PsiSqConj = conj(PsiSq);
	A11 = MCoef .* (PsiSqConj .* Psi{1}.^2 - PsiSq .* PsiConj{1}.^2);
	A12 = MCoef .* (PsiSqConj .* Psi{1} .* Psi{2} - PsiSq .* PsiConj{1} .* PsiConj{2});
	A13 = MCoef .* (PsiSqConj .* Psi{1} .* Psi{3} - PsiSq .* PsiConj{1} .* PsiConj{3});
	A22 = MCoef .* (PsiSqConj .* Psi{2}.^2 - PsiSq .* PsiConj{2}.^2);
	A23 = MCoef .* (PsiSqConj .* Psi{2} .* Psi{3} - PsiSq .* PsiConj{2} .* PsiConj{3});
	A33 = MCoef .* (PsiSqConj .* Psi{3}.^2 - PsiSq .* PsiConj{3}.^2);
	
	p1 = A11.^2 + A13.^2 + A23^2;
	q = A11 .* A22 .* A33;
	p2 = (A11 - q).^2 + (A22 - q).^2 + (A33 - q).^2 + 2*p1;
	p = sqrt(p2 / 6);
	B11 = (A11 - q) ./ p;
	B12 = A12 ./ p;
	B13 = A13 ./ p;
	B22 = (A22 - q) ./ p;
	B23 = A23 ./ p;
	B33 = (A33 - q) ./ p;

	r = 0.5 * (B11.*B22.*B33 + 2 * (B12.*B23.*B13) - (A13.^2).*B22 - (A12.^2).*B33 - (A23.^2).*B11);
	phi = acos(r) / 3;

	lambda1 = q + 2*p.*cos(phi);
	lambda3 = q + 2*p.*cos(phi + (2*pi/3));
	lambda2 = 3 * q - lambda1 - lambda3;

	expl1 = exp(lambda1);
	expl2 = exp(lambda2);
	expl3 = exp(lambda3);

	% evec1 = cross()

	
	% Rho = getRho(Psi);
	% RhoSq = Rho .^ 2;
	% PsiSq = dotVectorArray(Psi, Psi);
	% PsiSqAbsSq = abs(PsiSq) .^ 2;
	% Lambda1 = sqrt(PsiSqAbsSq.^2 - PsiSqAbsSq .* RhoSq);
	% Lambda2 = -Lambda1;
	% Lambda3 = 0;
	% gamma = (-Rho .* PsiSq) .* Lambda1 ./ (1 + PsiSqAbsSq .* Lambda1);
	% gamma = (-Rho .* PsiSq) .* Lambda1 ./ (1 + PsiSqAbsSq .* Lambda1);
end
function [o1, o2, o3] = updatePointPsi(psi1, psi2, psi3, dt, simConfig)
	% psiSq = psi1^2 + psi2^2 + psi3^2;
	psi = [psi1; psi2; psi3];
	psiSq = psi1^2 + psi2^2 + psi3^2;
	psiSqDag = conj(psiSq);
	psiDag = conj(psi);
	matr = (psi * transpose(psi) * psiSqDag - psiDag * transpose(psiDag) * psiSq) * (simConfig.siCoef^2 * dt^2 * -0.5) + (psiDag * transpose(psi)) * (simConfig.siCoef * dt * -1i);
	psiUpdated = expm(matr) * psi;
	o1 = psiUpdated(1);
	o2 = psiUpdated(2);
	o3 = psiUpdated(3);
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
	Rho = getRho(PsiForOp);
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