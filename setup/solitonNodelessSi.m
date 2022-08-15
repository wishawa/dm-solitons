% function Out = solitonNodelessSi(ctr, rc, epsilon, simConfig)
% 	Out = createRepeatingSolitons(simConfig, @(sp) createInner(sp, ctr, rc, epsilon, simConfig), 1);
% end

function Out = solitonNodelessSi(ctr, amp, epsilon, simConfig)
	epsilon = epsilon / norm(epsilon);
	R = getR(ctr, simConfig);
	Pro = amp * sqrt(0.998 / (1. + 0.377 * R.^2 / amp));
	Out = cell(1, 3);
	for j = 1:3
		Out{j} = Pro * epsilon(j);
	end
end
