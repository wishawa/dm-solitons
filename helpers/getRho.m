function Rho = getRho(Psi)
%myFun - Description
%
% Syntax: Rho = getRho(Psi)
%
% Long description
	
Rho = zeros(size(Psi{1}), 'double');
for j = 1:3
	psj = Psi{j};
	Rho = Rho + abs(psj).^2;
end
end