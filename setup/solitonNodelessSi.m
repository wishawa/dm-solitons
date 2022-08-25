% function Out = solitonNodelessSi(ctr, rc, epsilon, simConfig)
% 	Out = createRepeatingSolitons(simConfig, @(sp) createInner(sp, ctr, rc, epsilon, simConfig), 1);
% end

function Out = solitonNodelessSi(ctr, r95, epsilon, simConfig)
	amp = 20.9024 / r95^2;
	epsilon = epsilon / norm(epsilon);
	R = getR(ctr, simConfig);
	if simConfig.useNoSiProfile
		Profile = amp * 0.99915 / (1. + 0.037653 * amp * R.^2).^4;
	else
		polarization = 3 - norm(cross(1i * epsilon, conj(epsilon)));
		ProfileMain = 0.998309 / (1. + 0.037653 * amp * R.^2).^8;
		ProfileCorrection = 0.089017 * amp^2 * R.^2 * polarization * simConfig.lambda ./ (1. + 0.050255 * amp * R.^2).^8;
		Profile = amp * sqrt(ProfileMain + ProfileCorrection);
	end
	Out = cell(1, 3);
	for j = 1:3
		Out{j} = Profile * epsilon(j);
	end
end
