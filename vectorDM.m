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

% Derived Constants (don't change)
simConfig.m = simConfig.m22 * 8.96215327e-89;	% 10^-22 eV / c^2 / mass of sun
simConfig.m_per_hbar = simConfig.m / hbar;
simConfig.dx = simConfig.Lbox / simConfig.N;
simConfig.siCoef = simConfig.lambda / (4 * simConfig.m * c * simConfig.m_per_hbar^2);

% Display Parameters
simConfig.plotEvery = 8;
simConfig.plotGridBoxSize = 8;

% Simulation Parameters
simConfig.totalIterations = 12000;
simConfig.snapEvery = 400;

simConfig.plotGridBoxSize = 16;
% nSols = 8;
% simConfig.positions = rand(nSols, 3) .* [simConfig.N simConfig.N simConfig.N] - simConfig.N/2;
% simConfig.sizes = rand(nSols, 1) * 6 + 0.8;
% simConfig.spins = rand(nSols, 3) + rand(nSols, 3)*1i;

% simConfig.lambda = 1E-84;
% simulate("outputs/2022-07-19/8-solitons-random-144-repulsive", simConfig);
% simConfig.lambda = 0;
% simulate("outputs/2022-07-19/8-solitons-random-144-nosi", simConfig);
% simConfig.lambda = -1E-84;
% simulate("outputs/2022-07-19/8-solitons-random-144-attractive", simConfig);
simConfig.lambda = -1E-84;
simConfig.ctrs = [0 0 0];
simConfig.sizes = [0.75];
simConfig.epsilons = [1 1i 0];
simulate("outputs/2022-07-21/1-soliton-0.75,circular", simConfig);
simConfig.epsilons = [1 1 0];
simulate("outputs/2022-07-21/1-soliton-0.75,linear", simConfig);
% for i = 1:5
% 	pa = "outputs/2022-07-20/8-solitons-random-128-attractive-run-" + i;
% 	pn = "outputs/2022-07-21/8-solitons-random-128-nosi-run-" + i;
% 	simConfig = load(pa + "/simConfig.mat").simConfig;
% 	simConfig.lambda = 0E-84;
% 	simulate(pn, simConfig);
% end

function simulate(savename, simConfig)
	arguments
		savename string
		simConfig struct
	end

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

