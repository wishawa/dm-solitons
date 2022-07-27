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
addpath('setup/')		% for specifying spatial properties of the initial field
addpath('math/')

fftw('planner', 'measure');

% Constants
hbar = 1.71818131e-87;		% hbar / (mass of sun * (km/s) * kpc)
G = 4.3022682e-6;			% G/((km/s)^2*kpc/mass of sun)
c = 299792.458;				% c / (km/s)

simConfig = struct;

% Constants
simConfig.hbar = hbar;
simConfig.G = G;
simConfig.c = c;

% Chosen Constants
simConfig.m22 = 100;
simConfig.Lbox = 100.0;
simConfig.N = 128;
simConfig.lambda = -1E-84;

% Debug Parameters
simConfig.useSponge = false;
simConfig.dtOver = 1;
simConfig.doDrift = true;
simConfig.doScalarKick = true;
simConfig.doVectorKick = true;
simConfig.doVectorCorrection = true;

% Display Parameters
simConfig.plotEvery = 40;
simConfig.plotGridBoxSize = 16;

% Simulation Parameters
simConfig.totalIterations = 6400;
simConfig.snapEvery = 800;

for i = 1:5
	[simConfig.ctrs, simConfig.sizes, simConfig.epsilons] = randomSolitonsConfigs(8, 1.0, 4.0, simConfig.Lbox);
	simConfig.lambda = 1E-84;
	simulate("outputs/2022-07-27/8-solitons-random-96-repulsive-run-" + i, simConfig);
	simConfig.lambda = 0;
	simulate("outputs/2022-07-27/8-solitons-random-96-nosi-run-" + i, simConfig);
	simConfig.lambda = -1E-84;
	simulate("outputs/2022-07-27/8-solitons-random-96-attractive-run-" + i, simConfig);
end

function simulate(savename, simConfig)
	arguments
		savename string
		simConfig struct
	end

	simConfig.m = simConfig.m22 * 8.96215327e-89;
	simConfig.m_per_hbar = simConfig.m / simConfig.hbar;
	simConfig.dx = simConfig.Lbox / simConfig.N;
	simConfig.siCoef = simConfig.lambda / (4 * simConfig.m * simConfig.c * simConfig.m_per_hbar^2);

	mkdir(savename);

	iterations = simConfig.totalIterations;
	N = simConfig.N;
	Lbox = simConfig.Lbox;

	save(sprintf("%s/simConfig.mat", savename), 'simConfig');

	Psi = solitonsFromConfigs(simConfig);
	Rho = getRho(Psi, simConfig);
	totalMass = getTotalMass(Rho, simConfig);
	rhobar = totalMass / Lbox^3;

	klin = (-N/2:N/2-1)' * (2*pi/Lbox);
	[k1, k2, k3] = meshgrid(klin, klin, klin);
	kSq = fftshift(k1.^2 + k2.^2 + k3.^2);
	kSqNonzero = kSq + (kSq == 0);

	t = 0;
	i = 0;

	cflSchrodinger = (simConfig.m_per_hbar / 6) * simConfig.dx^2;

	displayer = SimulationDisplayer(simConfig, savename);
	displayer.displayStep(Psi, t);

	Rho = getRho(Psi, simConfig);
	VGrav = getGravPotential(Rho, rhobar, kSqNonzero, simConfig);
	VSiScalar = getSiScalarPotential(Rho, simConfig);
	VScalar = VGrav + VSiScalar;

	while i < iterations
		% Time Conditions
		cflNonlinear = pi / (max(abs(VScalar), [], 'all'));
		dt = min(cflSchrodinger, cflNonlinear) / simConfig.dtOver;

		% Drift
		if (simConfig.doDrift)
			Psi = stepDrift(Psi, kSq, dt / 2, simConfig);
		end

		% Update Potentials
		Rho = getRho(Psi, simConfig);
		VGrav = getGravPotential(Rho, rhobar, kSqNonzero, simConfig);
		VSiScalar = getSiScalarPotential(Rho, simConfig);
		VScalar = VGrav + VSiScalar;

		% Kick
		if (simConfig.doScalarKick)
			Psi = stepKickScalar(Psi, VScalar, dt/2);
		end
		if (simConfig.doVectorKick)
			Psi = stepKickVector(Psi, dt, simConfig);
		end
		if (simConfig.doScalarKick)
			Psi = stepKickScalar(Psi, VScalar, dt/2);
		end

		% Drift
		if (simConfig.doDrift)
			Psi = stepDrift(Psi, kSq, dt / 2, simConfig);
		end

		% Absorb
		if (simConfig.useSponge)
			centerRange = floor(N / 8):ceil(N * 7 / 8);
			for j = 1:3
				CenterPsi = Psi{j}(centerRange, centerRange, centerRange);
				Psi{j} = zeros(size(Psi{j}));
				Psi{j}(centerRange, centerRange, centerRange) = CenterPsi;
			end
		end

		t = t + dt;

		if rem(i, simConfig.snapEvery) == 0
			save(sprintf("%s/snap-Psi-%d-%.2f.mat", savename, i, t), 'Psi');
		end

		% Display
		displayer.displayStep(Psi, t);

		i = i + 1;
	end
	displayer.finish();
	save(sprintf("%s/snap-Psi-%d-%.2f.mat", savename, i, t), 'Psi');
end

