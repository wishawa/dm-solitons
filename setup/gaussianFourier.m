function Out1 = gaussianFourier(sigma, targetDensity, simConfig)
	N = simConfig.N;
	Lbox = simConfig.Lbox;
	klin = (-N/2:N/2-1)' * (2*pi/Lbox);
	[k1, k2, k3] = meshgrid(klin, klin, klin);
	kSq = fftshift(k1.^2 + k2.^2 + k3.^2);

	% FourierPsi1 = sqrt(5000 * exp(-kSq));
	% Out1 = ifftn(FourierPsi1) .* exp(1i * 2 * pi * rand(N, N, N));
	Out1 = exp(1i * 2 * pi * rand(size(kSq)));
	Out1 = Out1 .* sqrt(exp(-kSq / (2 * sigma^2)));
	Out1 = ifftn(Out1);
	Out1 = Out1 * sqrt(targetDensity / (sum(abs(Out1(:)).^2) * simConfig.dx^3 / simConfig.Lbox^3));
end