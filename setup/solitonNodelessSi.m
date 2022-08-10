% function Out = solitonNodelessSi(ctr, rc, epsilon, simConfig)
% 	Out = createRepeatingSolitons(simConfig, @(sp) createInner(sp, ctr, rc, epsilon, simConfig), 1);
% end

function Out = solitonNodelessSi(ctr, rc, epsilon, simConfig)
	% lambda = simConfig.lambda;
	m22 = simConfig.m22;
	epsilon = epsilon / norm(epsilon);
	% pol = 3.0 - norm(cross(1i * epsilon, conj(epsilon)));
	% rceff = rc * (1. + 3./8. * lambda * pol * 2.74E90 / m22^4 / rc^2);
	rceff = rc;
	rho0 = 1.9E7 * m22^-2 * rceff^-4;
	R = getR(ctr, simConfig);
	Dns = sqrt(rho0 ./ (1 + 0.091 * (R / rceff).^2).^8);
	Out = cell(1, 3);
	for j = 1:3
		Out{j} = Dns * epsilon(j);
	end
end
