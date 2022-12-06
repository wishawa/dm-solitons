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

simPath = "out_remote/2022-12-05/condensation,sigma=1,L=50,density=0.00016,lambda=-10";

% simulate("outputs/_testbed", 5.21, 43, simConfig);
simulate(simPath, 400000, 366210.93);

function simulate(simPath, iStart, tStart)

	simConfig = load(sprintf("%s/simConfig.mat", simPath)).simConfig;
	simConfig.totalIterations = 600000;

	rng(1234);

	iterations = simConfig.totalIterations;
	N = simConfig.N;
	Lbox = simConfig.Lbox;

	% Psi = solitonsFromConfigs(simConfig);
	% Psi = {gaussianFourier(sigma, targetDensity, simConfig), zeros(N, N, N), zeros(N, N, N)};
	Psi = load(sprintf("%s/snap-Psi-%d-%.2f.mat", simPath, iStart, tStart)).Psi;

	Rho = getRho(Psi);
	totalMass = getTotalMass(Rho, simConfig);
	rhobar = totalMass / Lbox^3;

	klin = (-N/2:N/2-1)' * (2*pi/Lbox);
	[k1, k2, k3] = meshgrid(klin, klin, klin);
	kSq = fftshift(k1.^2 + k2.^2 + k3.^2);
	kSqNonzero = kSq + (kSq == 0);

	t = tStart;
	i = iStart;

	% cflSchrodinger = 2./pi * simConfig.dx^2;
	% cflSchrodinger = 1./2. * simConfig.dx^2;
	% cflSchrodinger = 4 * simConfig.dx^2;
	cflSchrodinger = 1.5 * simConfig.dx^2;

	mkdir(simPath + "/resumed")
	displayer = SimulationDisplayer(simConfig, sprintf("%s/resumed", simPath));
	displayer.displayStep(Psi, t);

	Rho = getRho(Psi);
	RhoMoved = Rho - rhobar;
	VGrav = getGravPotential(RhoMoved, kSqNonzero);
	VSiScalar = getSiScalarPotential(RhoMoved, simConfig);
	VScalar = VGrav + VSiScalar;


	% save(sprintf("%s/snap-Psi-%d-%.2f.mat", savename, i, t), 'Psi');
	while i < iterations
        % tic;

		% Time Conditions
		cflNonlinear = 1. / (max(abs(VScalar), [], 'all'));
		dt = min(cflSchrodinger, cflNonlinear) / simConfig.dtOver;
		% dt = cflSchrodinger / simConfig.dtOver;

		Psi = stepKickScalar(Psi, VScalar, dt/2.);
		Psi = stepDrift(Psi, kSq, dt);

		Rho = getRho(Psi);
		RhoMoved = Rho - rhobar;
		VGrav = getGravPotential(RhoMoved, kSqNonzero);
		VSiScalar = getSiScalarPotential(RhoMoved, simConfig);
		VScalar = VGrav + VSiScalar;

		Psi = stepKickScalar(Psi, VScalar, dt/2.);

		% % Drift
		% if (simConfig.doDrift)
		% 	Psi = stepDrift(Psi, kSq, dt / 2);
		% end

		% % Update Potentials
		% Rho = getRho(Psi);
		% VGrav = getGravPotential(Rho, rhobar, kSqNonzero);
		% VSiScalar = getSiScalarPotential(Rho, simConfig);
		% VScalar = VGrav + VSiScalar;

		% % Kick
		% if (simConfig.doScalarKick)
		% 	Psi = stepKickScalar(Psi, VScalar, dt/2);
		% end
		% if (simConfig.doVectorKick)
		% 	Psi = stepKickVector(Psi, dt, simConfig);
		% end
		% if (simConfig.doScalarKick)
		% 	Psi = stepKickScalar(Psi, VScalar, dt/2);
		% end

		% % Drift
		% if (simConfig.doDrift)
		% 	Psi = stepDrift(Psi, kSq, dt / 2);
		% end

		% % Absorb
		% if (simConfig.useSponge)
		% 	centerRange = floor(N / 8):ceil(N * 7 / 8);
		% 	for j = 1:3
		% 		CenterPsi = Psi{j}(centerRange, centerRange, centerRange);
		% 		Psi{j} = zeros(size(Psi{j}));
		% 		Psi{j}(centerRange, centerRange, centerRange) = CenterPsi;
		% 	end
		% end

		t = t + dt;
		i = i + 1;

		if (rem(i, simConfig.snapEvery) == 0 || ((i > simConfig.totalIterations - simConfig.endSnapsIterations) && rem(i, simConfig.endSnapEvery) == 0))
			save(sprintf("%s/snap-Psi-%d-%.2f.mat", simPath, i, t), 'Psi');
		end

		% Display
		displayer.displayStep(Psi, t);

        % toc
	end
	displayer.finish();
	save(sprintf("%s/snap-Psi-%d-%.2f.mat", simPath, i, t), 'Psi');
end

