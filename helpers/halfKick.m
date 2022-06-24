function NewPsi = halfKick(Psi, VScalar, Rho, PsiForVsi, dt, simConsts)
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
N = simConsts.N;
for j = 1:3
	buffer = zeros(N, N, N);
	for k = 1:3
		delta = double(j == k);
		buffer = buffer + (delta + MCoef .* conj(PsiForVsi{j}) .* PsiForVsi{k}) .* BPsi{k};
	end
	NewPsi{j} = buffer;
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