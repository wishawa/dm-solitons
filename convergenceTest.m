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

% simConfig = struct;

% % Chosen Constants
% simConfig.Lbox = 400.0;
% simConfig.N = 96;
% simConfig.lambda = 0;

% % Debug Parameters
% simConfig.useSponge = false;
% simConfig.dtOver = 1;
% simConfig.doDrift = true;
% simConfig.useNoSiProfile = false;
% simConfig.doScalarKick = true;
% simConfig.doVectorKick = true;
% simConfig.doVectorCorrection = true;

% % Display Parameters
% simConfig.plotEvery = 8;
% simConfig.plotGridBoxSize = 16;

% % Simulation Parameters
% simConfig.totalIterations = 6000;
% simConfig.snapEvery = 100;
% simConfig.endSnapEvery = 100;
% simConfig.endSnapsIterations = 0;

% [simConfig.ctrs, simConfig.r95s, simConfig.epsilons] = randomSolitonsConfigs(5, 20.0, 40.0, simConfig.Lbox);
% simConfig.epsilons(1, :) = randomEpsilon(1);
% simConfig.epsilons(2, :) = randomEpsilon(0);
% for j = 1:2
% 	% for i = [1, 2, 4, 8]
% 	% 	simConfig.dtOver = i;
% 	% 	simConfig.totalIterations = 6000 * i;
% 	% 	simConfig.snapEvery = 100 * i;
% 	% 	simConfig.plotEvery = 8 * i;
% 	simulate(simConfig);
% 	% end
% end
% simConfig.ctrs = [0 0 0];
% simConfig.r95s = [20.0];
% simConfig.epsilons = [1 1i 0];
% simulate("outputs/_testbed", simConfig);
simConfig = load("outputs/simConfig.mat").simConfig;
simulate(simConfig);

function simulate(simConfig)
	arguments
		simConfig struct
	end
	simConfig.dtOver = 1;
	stepper1 = SimulationStepper(simConfig);
	simConfig.dtOver = 2;
	stepper2 = SimulationStepper(simConfig);
	simConfig.dtOver = 4;
	stepper4 = SimulationStepper(simConfig);


	% obj.displayer = SimulationDisplayer(simConfig, obj.saveName);
	% obj.displayer.displayStep(Psi, saveName);

	% save(sprintf("%s/snap-Psi-%d-%.2f.mat", saveName, obj.iter, obj.time), 'Psi');
	for i = 1:100

		stepper1.step();

		stepper2.step();
		stepper2.step();

		stepper4.step();
		stepper4.step();
		stepper4.step();
		stepper4.step();

		convergenceFactor = log2(convergenceDifference(stepper1.PsiO, stepper2.PsiO) ./ convergenceDifference(stepper2.PsiO, stepper4.PsiO));
		disp(prctile(convergenceFactor(:), [5, 25, 50, 75, 95]));

		% if (rem(obj.iter, simConfig.snapEvery) == 0 || ((obj.iter > simConfig.totalIterations - simConfig.endSnapsIterations) && rem(obj.iter, simConfig.endSnapEvery) == 0))
		% 	save(sprintf("%s/snap-Psi-%d-%.2f.mat", obj.saveName, obj.iter, obj.time), 'Psi');
		% end

		% % Display
		% obj.displayer.displayStep(Psi, t);
	end

	% obj.displayer.finish();
	% save(sprintf("%s/snap-Psi-%d-%.2f.mat", obj.saveName, obj.iter, obj.time), 'Psi');
end
function Dif = convergenceDifference(PsiL, PsiR)
	Dif = sqrt(getRho(cellfun(@(l, r) l - r, PsiL, PsiR, 'UniformOutput', false)));
end

