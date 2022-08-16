% function Out = solitonNodelessSi(ctr, rc, epsilon, simConfig)
% 	Out = createRepeatingSolitons(simConfig, @(sp) createInner(sp, ctr, rc, epsilon, simConfig), 1);
% end

function Out = solitonNodelessSi(ctr, r95, epsilon, simConfig)
	amp = 20.9024 / r95^2;
	epsilon = epsilon / norm(epsilon);
	R = getR(ctr, simConfig);
	Pro = amp * 0.99915 / (1. + 0.037653 * amp * R.^2).^4;
	Out = cell(1, 3);
	for j = 1:3
		Out{j} = Pro * epsilon(j);
	end
end
