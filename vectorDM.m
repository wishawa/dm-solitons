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
simConfig.endSnapsIterations = 800;

simConfig.endSnapEvery = 100;
simConfig.endSnapsIterations = 800;

% for i = 1:3:24
% 	[simConfig.ctrs, simConfig.r95s, simConfig.epsilons] = randomSolitonsConfigs(8, 2.0, 4.0, simConfig.Lbox);


% 	% simConfig = load("out_remote/2022-07-30/8-solitons-random-128-repulsive-run-" + i + "/simConfig.mat").simConfig;
% 	simConfig.lambda = 2E-84;
% 	simulate("out_remote/2022-08-15/8-solitons-random-128-strong-repulsive-run-" + i, simConfig);

% 	simConfig.lambda = 1E-84;
% 	simulate("out_remote/2022-08-15/8-solitons-random-128-repulsive-run-" + i, simConfig);

% 	simConfig.lambda = 5E-85;
% 	simulate("out_remote/2022-08-15/8-solitons-random-128-weak-repulsive-run-" + i, simConfig);

% 	% simConfig = load("out_remote/2022-07-30/8-solitons-random-128-nosi-run-" + i + "/simConfig.mat").simConfig;
% 	simConfig.lambda = 0;
% 	simulate("out_remote/2022-08-15/8-solitons-random-128-nosi-run-" + i, simConfig);

% 	simConfig.lambda = -5E-85;
% 	simulate("out_remote/2022-08-15/8-solitons-random-128-weak-attractive-run-" + i, simConfig);

% 	% simConfig = load("out_remote/2022-07-30/8-solitons-random-128-attractive-run-" + i + "/simConfig.mat").simConfig;
% 	simConfig.lambda = -1E-84;
% 	simulate("out_remote/2022-08-15/8-solitons-random-128-attractive-run-" + i, simConfig);

% 	simConfig.lambda = -2E-84;
% 	simulate("out_remote/2022-08-15/8-solitons-random-128-strong-attractive-run-" + i, simConfig);
% end
simConfig.lambda = 0;
simConfig.N = 128;
simConfig.plotEvery = 1;
simConfig.totalIterations = 8000;
simConfig.ctrs = [0 0 0];
simConfig.r95s = [25];
simConfig.epsilons = [1 0 0];
simConfig.doVectorCorrection = false;
simulate("outputs/test_unitless", simConfig);

function simulate(savename, simConfig)
	arguments
		savename string
		simConfig struct
	end

	simConfig.dx = simConfig.Lbox / simConfig.N;

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

	cflSchrodinger = 1./6. * simConfig.dx^2;

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

