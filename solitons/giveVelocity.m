function Out = giveVelocity(Spaces, Soliton, Velocity, simConsts)
%myFun - Description
%
% Syntax: Out = giveVelocity(soliton)
%
% Long description
	Out = cell(1, 3);
	for j = 1:3
		Out{j} = Soliton{j};
		for k = 1:3
			Out{j} = Out{j} .* exp(1i * simConsts.m_per_hbar * Spaces{k} .* Velocity(k));
		end
	end
end