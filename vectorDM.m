%% Rohith Karur (2022), UC Berkeley, Philip Mocz (2021), Princeton University
% Merge solitons with Vector Dark Matter Formulation

%% Code modified for 2 soliton purpose by Rohith Karur 2021
% The purpose of this code is to simulate 2 soliton interactions

% Internal units:
% [L] = kpc
% [M] = Msun
% [E] = Msun (km/s)^2

% scaling
% {x, t, rho, m} --> {a x, b t, b^-2 rho, a^-2 b m}

% v = (hbar / m) * grad(phase)


addpath('helpers/')			% functions for extracting energies, potential etc. 
addpath('solitons/')		% for specifying spatial properties of the initial field

fftw('planner', 'measure');

simulate(100, 100.0, 96, 1E-84, @(Spaces, simConsts) {...
	giveVelocity(Spaces, solitonNodelessSi(Spaces, [0 -20 0], 4.0, [1 1i 0], simConsts), [0 0.04 0], simConsts),...
	giveVelocity(Spaces, solitonNodelessSi(Spaces, [0 20 0], 4.0, [1 0 0], simConsts), [0 -0.04 0], simConsts),...
}, 8, 12, "outputs/2022-07-06/collide-repulsive-linear-circular", 10000, false);


function simulate(m22, Lbox, N, lambda, createSolitons, snapEvery, gridEvery, savename, iterations, useSponge)
	arguments
		m22 double
		Lbox double
		N double
		lambda double
		createSolitons
		snapEvery int32
		gridEvery int32
		savename string
		iterations int32
		useSponge logical
	end

	mkdir(savename);

	% Constants
	hbar = 1.71818131e-87;		% hbar / (mass of sun * (km/s) * kpc)
	G = 4.3022682e-6;			% G/((km/s)^2*kpc/mass of sun)
	c = 299792.458;				% c / (km/s)

	% Simulation Constants
	simConsts = struct;
	simConsts.totalIterations = iterations;

	simConsts.hbar = hbar;
	simConsts.G = G;
	simConsts.c = c;

	simConsts.m22 = m22;
	simConsts.Lbox = Lbox;
	simConsts.N = N;
	simConsts.lambda = lambda;

	m = m22 * 8.96215327e-89;	% 10^-22 eV / c^2 / mass of sun
	m_per_hbar = m / hbar;
	dx = Lbox / N;
	siCoef = lambda / (4 * m * c * m_per_hbar^2);

	simConsts.m = m;
	simConsts.m_per_hbar = m_per_hbar;
	simConsts.dx = dx;
	simConsts.siCoef = siCoef;

	Psi = createRepeatingSolitons(simConsts, createSolitons);

	Rho = getRho(Psi, simConsts);
	totalMass = getTotalMass(Rho, simConsts);
	rhobar = totalMass / Lbox^3;

	klin = (-N/2:N/2-1)' * (2*pi/Lbox);
	[k1, k2, k3] = meshgrid(klin, klin, klin);
	kSq = fftshift(k1.^2 + k2.^2 + k3.^2);
	kSqNonzero = kSq + (kSq == 0);

	t = 0;
	i = 0;

	cflSchrodinger = (m_per_hbar / 6) * dx^2;

	% tic;
	% lastToc = toc;

	displayer = SimulationDisplayer(simConsts, savename, snapEvery, gridEvery);
	displayer.displayStep(Psi, t);

	Rho = getRho(Psi, simConsts);
	VGrav = getGravPotential(Rho, rhobar, kSqNonzero, simConsts);
	VSiScalar = getSiScalarPotential(Rho, simConsts);
	VScalar = VGrav + VSiScalar;

	while i < iterations
		% Time Conditions
		cflNonlinear = pi / (max(abs(VScalar), [], 'all'));
		dt = min(cflSchrodinger, cflNonlinear);

		% Drift
		Psi = stepDrift(Psi, kSq, dt / 2, simConsts);

		% Update Potentials
		Rho = getRho(Psi, simConsts);
		VGrav = getGravPotential(Rho, rhobar, kSqNonzero, simConsts);
		VSiScalar = getSiScalarPotential(Rho, simConsts);
		VScalar = VGrav + VSiScalar;

		% Kick
		Psi = stepKickScalar(Psi, VScalar, dt/2);
		Psi = stepKickVector(Psi, Psi, Rho, dt, simConsts);
		Psi = stepKickScalar(Psi, VScalar, dt/2);

		% Drift
		Psi = stepDrift(Psi, kSq, dt / 2, simConsts);

		% Absorb
		if (useSponge)
			centerRange = floor(N / 8):ceil(N * 7 / 8);
			for j = 1:3
				CenterPsi = Psi{j}(centerRange, centerRange, centerRange);
				Psi{j} = zeros(size(Psi{j}));
				Psi{j}(centerRange, centerRange, centerRange) = CenterPsi;
			end
		end

		t = t + dt;

		if rem(i, 250) == 0
			save(sprintf("%s/snap-Psi-%d-%.2f.mat", savename, i, t), 'Psi');
		end

		% Display
		displayer.displayStep(Psi, t);

		% tocDif = toc - lastToc;
		% lastToc = toc;
		% disp("time taken: " + tocDif);
		i = i + 1;
	end
	displayer.finish();
end