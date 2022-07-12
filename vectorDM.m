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

% Constants
hbar = 1.71818131e-87;		% hbar / (mass of sun * (km/s) * kpc)
G = 4.3022682e-6;			% G/((km/s)^2*kpc/mass of sun)
c = 299792.458;				% c / (km/s)

simConsts = struct;

% Constants
simConsts.hbar = hbar;
simConsts.G = G;
simConsts.c = c;

% Chosen Constants
simConsts.m22 = 100;
simConsts.Lbox = 100.0;
simConsts.N = 96;
simConsts.lambda = -1E-84;

% Debug Parameters
simConsts.useSponge = false;
simConsts.dtOver = i;
simConsts.doDrift = true;
simConsts.doScalarKick = true;
simConsts.doVectorKick = true;

% Derived Constants
simConsts.m = m22 * 8.96215327e-89;	% 10^-22 eV / c^2 / mass of sun
simConsts.m_per_hbar = simConsts.m / hbar;
simConsts.dx = simConsts.Lbox / simConsts.N;
simConsts.siCoef = simConsts.lambda / (4 * simConsts.m * c * simConsts.m_per_hbar^2);

% Display Parameters
simConsts.snapEvery = 8;
simConsts.gridResolution = 8;

% Simulation Parameters
simConsts.totalIterations = 2000;

for i = 1:4
	simConsts.totalIterations = 4000 * i;
	simConsts.doVectorKick = true;
	Psi = giveVelocity(solitonNodelessSi([0 0 0], 4.0, [-1 1i 0], simConsts), [0.5 0 0], simConsts);
	simulate(simConsts, Psi, "outputs/2022-07-13/single-fast-moving-soliton-96,long,dto" + i);
end
for i = 1:4
	simConsts.totalIterations = 4000 * i;
	simConsts.doVectorKick = false;
	Psi = giveVelocity(solitonNodelessSi([0 0 0], 4.0, [-1 1i 0], simConsts), [0.5 0 0], simConsts);
	simulate(simConsts, Psi, "outputs/2022-07-13/single-fast-moving-soliton-96,long,novsi,dto" + i);
end

function simulate(simConsts, Psi, savename)
	arguments
		simConsts struct
		Psi
		savename string
	end

	mkdir(savename);

	iterations = simConsts.totalIterations;
	N = simConsts.N;
	Lbox = simConsts.Lbox;

	save(sprintf("%s/simConsts.mat", savename), 'simConsts');

	Rho = getRho(Psi, simConsts);
	totalMass = getTotalMass(Rho, simConsts);
	rhobar = totalMass / Lbox^3;

	klin = (-N/2:N/2-1)' * (2*pi/Lbox);
	[k1, k2, k3] = meshgrid(klin, klin, klin);
	kSq = fftshift(k1.^2 + k2.^2 + k3.^2);
	kSqNonzero = kSq + (kSq == 0);

	t = 0;
	i = 0;

	cflSchrodinger = (simConsts.m_per_hbar / 6) * simConsts.dx^2;


	displayer = SimulationDisplayer(simConsts, savename);
	displayer.displayStep(Psi, t);

	Rho = getRho(Psi, simConsts);
	VGrav = getGravPotential(Rho, rhobar, kSqNonzero, simConsts);
	VSiScalar = getSiScalarPotential(Rho, simConsts);
	VScalar = VGrav + VSiScalar;

	while i < iterations
		% Time Conditions
		cflNonlinear = pi / (max(abs(VScalar), [], 'all'));
		dt = min(cflSchrodinger, cflNonlinear) / simConsts.dtOver;

		% Drift
		if (simConsts.doDrift)
			Psi = stepDrift(Psi, kSq, dt / 2, simConsts);
		end

		% Update Potentials
		Rho = getRho(Psi, simConsts);
		VGrav = getGravPotential(Rho, rhobar, kSqNonzero, simConsts);
		VSiScalar = getSiScalarPotential(Rho, simConsts);
		VScalar = VGrav + VSiScalar;

		% Kick
		if (simConsts.doScalarKick)
			Psi = stepKickScalar(Psi, VScalar, dt/2);
		end
		if (simConsts.doVectorKick)
			Psi = stepKickVector(Psi, Psi, Rho, dt, simConsts);
		end
		if (simConsts.doScalarKick)
			Psi = stepKickScalar(Psi, VScalar, dt/2);
		end

		% Drift
		if (simConsts.doDrift)
			Psi = stepDrift(Psi, kSq, dt / 2, simConsts);
		end

		% Absorb
		if (simConsts.useSponge)
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

		i = i + 1;
	end
	displayer.finish();
end