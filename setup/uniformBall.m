function Out = uniformBall(Spaces, m22, ctr, rc, dVec, simConfig)
	dVec = dVec / norm(dVec);
	R = zeros(size(Spaces{1}));
	for j = 1:3
		R = R + (Spaces{j} - ctr(j)).^2;
	end
	R = sqrt(R);
	rho0 = 2.0E8 * m22^-2 * rc^-4;
	Dns = sqrt(rho0 .* double(R < rc));
	Out = cell(1, 3);
	for j = 1:3
		Out{j} = Dns * dVec(j);
	end
end