function totalSpins = getTotalSpins(Spins)
%myFun - Description
%
% Syntax: totalSpins = getTotalSpin(Spins)
%
% Long description
totalSpins = cell(1, 3);
for j = 1:3
	totalSpins{j} = sum(Spins{j}, 'all');
end
end