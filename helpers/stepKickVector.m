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
	% 	matr = (psi * transpose(psi) * psiSqDag - psiDag * transpose(psiDag) * psiSq) * (simConfig.siCoef^2 * dt^2 * 0.5) + (psiDag * transpose(psi)) * (simConfig.siCoef * dt * -1i);
	% 	psiUpdated = expm(matr) * psi;
	% 	Psi1(l) = psiUpdated(1);
	% 	Psi2(l) = psiUpdated(2);
	% 	Psi3(l) = psiUpdated(3);
	% end
	% Psi = {Psi1, Psi2, Psi3};

	OrigPsi = Psi;
	Psi = kickMainTerm(Psi, OrigPsi, dt / 2, simConfig);
	% if (simConfig.doVectorCorrection)
	% 	OrigPsiConj = {conj(OrigPsi{1}), conj(OrigPsi{2}), conj(OrigPsi{3})};
	% 	Psi = kickCorrectionTerm(Psi, OrigPsi, 1, dt, simConfig);
	% 	Psi = kickCorrectionTerm(Psi, OrigPsiConj, -1, dt, simConfig);
	% end
	if (simConfig.doVectorCorrection)
		Psi = kickCorrectionNew(Psi, OrigPsi, dt, simConfig);
	end

	Psi = kickMainTerm(Psi, OrigPsi, dt / 2, simConfig);

	% if (simConfig.doVectorCorrection)
	% 	Psi = kickCorrectionTerm(Psi, OrigPsiConj, -1, dt / 2, simConfig);
	% 	Psi = kickCorrectionTerm(Psi, OrigPsi, 1, dt / 2, simConfig);
	% end
end
function EigVec = evecFromLambda(A11, A12, A13, A22, A23, A33, lambda)
	M11 = A11 - lambda;
	M22 = A22 - lambda;
	M33 = A33 - lambda;
	R1 = [M11, A12, A13];
	R2 = [A12, M22, A23];
	R3 = [A13, A23, M33];
	nLattice = length(lambda);
	Cro = zeros(nLattice, 3, 3, 'double');
	Cro(:, :, 1) = cross(R2, R3);
	Cro(:, :, 2) = cross(R1, R3);
	Cro(:, :, 3) = cross(R1, R2);
	NormSq = zeros(nLattice, 3, 'double');
	CroSq = abs(Cro).^2;
	NormSq(:, 1) = sum(CroSq(:, :, 1), 2);
	NormSq(:, 2) = sum(CroSq(:, :, 2), 2);
	NormSq(:, 3) = sum(CroSq(:, :, 3), 2);
	[NormSqSum, Which] = max(NormSq, [], 2, 'linear');
	Co1 = Cro(:, 1, :);
	Co2 = Cro(:, 2, :);
	Co3 = Cro(:, 3, :);
	EigVec = [Co1(Which), Co2(Which), Co3(Which)];
	% disp(max(abs(EigVec - Cro(:, :, 1)), [], 'all') + " || " + max(abs(EigVec - Cro(:, :, 2)), [], 'all') + " || " + max(abs(EigVec - Cro(:, :, 1)), [], 'all'));
	% assert(max(abs(EigVec - Cro(:, :, 1)), [], 'all') < 1E-16 || max(abs(EigVec - Cro(:, :, 2)), [], 'all') < 1E-16 || max(abs(EigVec - Cro(:, :, 1)), [], 'all') < 1E-16);
	EigVec = EigVec ./ sqrt(NormSqSum);
	% EigVec = cross(R1, R2);
end
function NewPsi = kickCorrectionNew(TargetPsi, Psi, dt, simConfig)

	N = simConfig.N;

	PsiSq = dotVectorArray(Psi, Psi);
	PsiSq = PsiSq(:);
	PsiSqConj = conj(PsiSq);
	Psi = {Psi{1}(:), Psi{2}(:), Psi{3}(:)};
	PsiConj = {conj(Psi{1}), conj(Psi{2}), conj(Psi{3})};
	Rho = getRho(Psi);
	siCoef = simConfig.siCoef;
	ccoef = 0.5 * siCoef^2 * dt^2;

	A11 = real(1i * (PsiSq .* PsiConj{1}.^2				- PsiSqConj .* Psi{1}.^2 		));
	A12 = real(1i * (PsiSq .* PsiConj{1} .* PsiConj{2}	- PsiSqConj .* Psi{1} .* Psi{2} ));
	A13 = real(1i * (PsiSq .* PsiConj{1} .* PsiConj{3}	- PsiSqConj .* Psi{1} .* Psi{3} ));
	A22 = real(1i * (PsiSq .* PsiConj{2}.^2				- PsiSqConj .* Psi{2}.^2 		));
	A23 = real(1i * (PsiSq .* PsiConj{2} .* PsiConj{3}	- PsiSqConj .* Psi{2} .* Psi{3} ));
	A33 = real(1i * (PsiSq .* PsiConj{3}.^2				- PsiSqConj .* Psi{3}.^2 		));
	
	% tic;
	% p1 = A11.^2 + A13.^2 + A23.^2;
	% q = (A11 + A22 + A33) / 3.;
	% p2 = (A11 - q).^2 + (A22 - q).^2 + (A33 - q).^2 + 2*p1;
	% p = sqrt(p2 / 6);
	% B11 = (A11 - q) ./ p;
	% B12 = A12 ./ p;
	% B13 = A13 ./ p;
	% B22 = (A22 - q) ./ p;
	% B23 = A23 ./ p;
	% B33 = (A33 - q) ./ p;

	% r = 0.5 * (B11.*B22.*B33 + 2 * (B12.*B23.*B13) - (B13.^2).*B22 - (B12.^2).*B33 - (B23.^2).*B11);
	% r = min(max(r, -1.), 1.);
	% phi = acos(r) / 3;

	% oldlambda1 = q + 2*p.*cos(phi);
	% oldlambda3 = q + 2*p.*cos(phi + (2*pi/3));
	% oldlambda2 = 3 * q - oldlambda1 - oldlambda3;
	% toc
	lambda1 = 1i * abs(PsiSq) .* sqrt(abs(PsiSq).^2 - (Rho).^2);
	lambda2 = zeros(N^3, 1);
	lambda3 = -lambda1;
	% disp(max(abs(oldlambda1 - lambda1), [], 'all'));
	% disp(max(abs(oldlambda2 - lambda2), [], 'all'));
	% disp(max(abs(oldlambda3 - lambda3), [], 'all'));

	eD11 = exp(1i * ccoef * lambda1);
	eD22 = exp(1i * ccoef * lambda2);
	eD33 = exp(1i * ccoef * lambda3);

	Sx1 = evecFromLambda(A11, A12, A13, A22, A23, A33, lambda1);
	Sx2 = evecFromLambda(A11, A12, A13, A22, A23, A33, lambda2);
	Sx3 = evecFromLambda(A11, A12, A13, A22, A23, A33, lambda3);

	% Sx = {Sx1, Sx2, Sx3};
	% lambda = {lambda1, lambda2, lambda3};
	% for k = 1:3
	% 	assert(max(abs((A11 - lambda{k}) .* Sx{k}(:, 1) + A12 .* Sx{k}(:, 2) + A13 .* Sx{k}(:, 3)), [], 'all') < 1E-16);
	% 	assert(max(abs(A12 .* Sx{k}(:, 1) + (A22 - lambda{k}) .* Sx{k}(:, 2) + A23 .* Sx{k}(:, 3)), [], 'all') < 1E-16);
	% 	assert(max(abs(A13 .* Sx{k}(:, 1) + A23 .* Sx{k}(:, 2) + (A33 - lambda{k}) .* Sx{k}(:, 3)), [], 'all') < 1E-16);
	% end
	% midInd = sub2ind([N, N, N], N/2,N/2,N/2);
	% disp(TargetPsi{1}(midInd));

	TargetPsi = [TargetPsi{1}(:), TargetPsi{2}(:), TargetPsi{3}(:)];
	NewPsi = [dot(Sx1, TargetPsi, 2) .* eD11, dot(Sx2, TargetPsi, 2) .* eD22, dot(Sx3, TargetPsi, 2) .* eD33];
	NewPsi = NewPsi(:, 1).*Sx1 + NewPsi(:, 2).*Sx2 + NewPsi(:, 3).*Sx3;
	NewPsi = {...
		reshape(NewPsi(:, 1), [N, N, N]),...
		reshape(NewPsi(:, 2), [N, N, N]),...
		reshape(NewPsi(:, 3), [N, N, N]),...
	};
	
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
% function [o1, o2, o3] = updatePointPsi(psi1, psi2, psi3, dt, simConfig)
% 	% psiSq = psi1^2 + psi2^2 + psi3^2;
% 	psi = [psi1; psi2; psi3];
% 	psiSq = psi1^2 + psi2^2 + psi3^2;
% 	psiSqDag = conj(psiSq);
% 	psiDag = conj(psi);
% 	matr = (psi * transpose(psi) * psiSqDag - psiDag * transpose(psiDag) * psiSq) * (simConfig.siCoef^2 * dt^2 * 0.5) + (psiDag * transpose(psi)) * (simConfig.siCoef * dt * -1i);
% 	psiUpdated = expm(matr) * psi;
% 	o1 = psiUpdated(1);
% 	o2 = psiUpdated(2);
% 	o3 = psiUpdated(3);
% end
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