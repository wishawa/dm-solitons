function Rho = getRho(Psi, simConsts)
%myFun - Description
%
% Syntax: Rho = getRho(Psi)
%
% Long description
	
N = simConsts.N;
Rho = zeros(N, N, N, 'double');
for j = 1:3
	psj = Psi{j};
	Rho = Rho + abs(psj).^2;
end
end