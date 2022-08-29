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

simConfig = struct;

% Chosen Constants
simConfig.Lbox = 900.0;
simConfig.N = 128;
simConfig.lambda = -1;

% Debug Parameters
simConfig.useSponge = false;
simConfig.dtOver = 1;
simConfig.doDrift = true;
simConfig.doScalarKick = true;
simConfig.doVectorKick = true;
simConfig.doVectorCorrection = true;

% Display Parameters
simConfig.plotEvery = 20;
simConfig.plotGridBoxSize = 16;

% Simulation Parameters
simConfig.totalIterations = 12000;
simConfig.snapEvery = 4000;
simConfig.endSnapEvery = 100;
simConfig.endSnapsIterations = 0;

simConfig.lambda = 1.;
simConfig.N = 128;
simConfig.Lbox = 96.0;
simConfig.plotEvery = 8;
simConfig.totalIterations = 400;
simConfig.ctrs = [0. 0. 0.];
% simConfig.r95s = [16.];
simConfig.epsilons = [1. 1. 1.];
simConfig.solIdxs = [25];
% simConfig.epsilons = randomEpsilon()
simConfig.useNoSiProfile = false;
simConfig.doVectorKick = true;
simConfig.doVectorCorrection = true;

s = 0.0;
simConfig.useSponge = true;
% simConfig.epsilons = [1., exp(1i * asin(s)), 0.];
% simConfig.doVectorCorrection = false;
% simulate("out_remote/2022-08-29/1-soliton,repulsive,ind=25,spin=" + s, simConfig);
% simConfig.doVectorCorrection = true;
for s = 0.75:0.25:1.00
	simConfig.epsilons = [1., exp(1i * asin(s)), 0.];
	simulate("out_remote/2022-08-29/1-soliton,repulsive,coreplot,ind=25,spin=" + s, simConfig);
end

function simulate(savename, simConfig)
	arguments
		savename string
		simConfig struct
	end

	simConfig.dx = simConfig.Lbox / simConfig.N;
	simConfig.doVectorKick = simConfig.doVectorKick && simConfig.lambda ~= 0;

	mkdir(savename);

	iterations = simConfig.totalIterations;
	N = simConfig.N;
	Lbox = simConfig.Lbox;

	save(sprintf("%s/simConfig.mat", savename), 'simConfig');

	Psi = solitonsFromConfigs(simConfig);

	Rho = getRho(Psi);
	totalMass = getTotalMass(Rho, simConfig);
	rhobar = totalMass / Lbox^3;

	klin = (-N/2:N/2-1)' * (2*pi/Lbox);
	[k1, k2, k3] = meshgrid(klin, klin, klin);
	kSq = fftshift(k1.^2 + k2.^2 + k3.^2);
	kSqNonzero = kSq + (kSq == 0);

	t = 0;
	i = 0;

	cflSchrodinger = 2./pi * simConfig.dx^2;

	displayer = SimulationDisplayer(simConfig, savename);
	displayer.displayStep(Psi, t);

	Rho = getRho(Psi);
	VGrav = getGravPotential(Rho, rhobar, kSqNonzero);
	VSiScalar = getSiScalarPotential(Rho, simConfig);
	VScalar = VGrav + VSiScalar;
	while i < iterations
        tic;
		% Time Conditions
		cflNonlinear = pi / (max(abs(VScalar), [], 'all'));
		dt = min(cflSchrodinger, cflNonlinear) / simConfig.dtOver;

		% Drift
		if (simConfig.doDrift)
			Psi = stepDrift(Psi, kSq, dt / 2);
		end

		% Update Potentials
		Rho = getRho(Psi);
		VGrav = getGravPotential(Rho, rhobar, kSqNonzero);
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
			Psi = stepDrift(Psi, kSq, dt / 2);
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

		if (rem(i, simConfig.snapEvery) == 0 || ((i > simConfig.totalIterations - simConfig.endSnapsIterations) && rem(i, simConfig.endSnapEvery) == 0))
			save(sprintf("%s/snap-Psi-%d-%.2f.mat", savename, i, t), 'Psi');
		end

		% Display
		displayer.displayStep(Psi, t);

		i = i + 1;
        toc
	end
	displayer.finish();
	save(sprintf("%s/snap-Psi-%d-%.2f.mat", savename, i, t), 'Psi');
end

