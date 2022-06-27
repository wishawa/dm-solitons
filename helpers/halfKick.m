function NewPsi = halfKick(Psi, VScalar, PsiForVsi, Rho, dt, simConsts)
%myFun - Description
%
% Syntax: Psi = halfKick(Psi, Vgrav, VsiScalar, dt, simConsts)
%
% Long description

NewPsi = Psi;

for j = 1:3
	NewPsi{j} = exp(-1i * (dt / 4) * VScalar) .* NewPsi{j};
end

BPsi = NewPsi;

rhoMul = -1i * (dt / 2) * simConsts.siCoef;
MCoef = (exp(Rho * rhoMul) - 1) ./ Rho;
for j = 1:3
	for k = 1:3
		NewPsi{j} = NewPsi{j} + (MCoef .* conj(PsiForVsi{j}) .* PsiForVsi{k}) .* BPsi{k};
	end
end

% M = zeros(3, 3);
% for j = 1:3
% 	for k = 1:3
% 		M(j, k) = simConsts.siCoef * PsiForVsi{j}
% 	end
% end

for j = 1:3
	NewPsi{j} = exp(-1i * (dt / 4) * VScalar) .* NewPsi{j};
end
end