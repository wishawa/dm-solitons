function Out3 = solitonExact(ctr, solIdx, epsilon, simConfig)
	epsilon = epsilon / norm(epsilon);
	if (simConfig.lambda == 0)
		polarization = 0;
	else
		polarization = 3 - norm(cross(1i * epsilon, conj(epsilon)));
	end
	pr = load("profiles/" + simConfig.lambda + "_" + polarization + "_" + solIdx + ".mat");
	amplitude = pr.amplitude;
	normalizedCurve = pr.curve;
	curveDx = pr.dx / sqrt(amplitude);
	normalizedCurve = max(normalizedCurve, 0.);
	fullCurve = normalizedCurve * amplitude;
	R = getR(ctr, simConfig);
	RIdx = max(min(round(R(:) / curveDx), length(normalizedCurve)), 1);
	Profile = reshape(fullCurve(RIdx), size(R));
	Out3 = cell(1, 3);
	for j = 1:3
		Out3{j} = Profile * epsilon(j);
	end
end