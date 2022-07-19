function Out = solitonNodelessSi(ctr, rc, dVec, simConfig)
	Out = createRepeatingSolitons(simConfig, @(sp) createInner(sp, ctr, rc, dVec, simConfig), 1);
end

function Out = createInner(Spaces, ctr, rc, dVec, simConfig)
	lambda = simConfig.lambda;
	m22 = simConfig.m22;
	dVec = dVec / norm(dVec);
	if (simConfig.oldPol)
		pol = 1.0 + norm(cross(1i * dVec, conj(dVec)));
	else
		pol = 3.0 - norm(cross(1i * dVec, conj(dVec)));
	end
	rceff = rc * (1. + 3./8. * lambda * pol * 2.74E90 / m22^4 / rc^2);
	rho0 = 2.0E7 * m22^-2 * rceff^-4;
	R = getR(Spaces, ctr);
	Dns = sqrt(rho0 ./ (1 + 0.091 * (R / rceff).^2).^8);
	Out = cell(1, 3);
	for j = 1:3
		Out{j} = Dns * dVec(j);
	end
end