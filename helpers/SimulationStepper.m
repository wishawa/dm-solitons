classdef SimulationStepper < handle
	properties
		PsiO
		simConfigO
		avgRho
		kSq
		kSqNz
		displayer
		time
		iter
		cflSchrodinger
	end
	methods
		function obj = SimulationStepper(simConfig)
			simConfig.dx = simConfig.Lbox / simConfig.N;
			simConfig.doVectorKick = simConfig.doVectorKick && simConfig.lambda ~= 0;
			obj.simConfigO = simConfig;

			Psi = solitonsFromConfigs(simConfig);
			obj.PsiO = Psi;

			Rho = getRho(Psi);
			totalMass = getTotalMass(Rho, simConfig);
			obj.avgRho = totalMass / simConfig.Lbox^3;

			klin = (-simConfig.N/2:simConfig.N/2-1)' * (2*pi/simConfig.Lbox);
			[k1, k2, k3] = meshgrid(klin, klin, klin);
			obj.kSq = fftshift(k1.^2 + k2.^2 + k3.^2);
			obj.kSqNz = obj.kSq + (obj.kSq == 0);

			obj.time = 0;
			obj.iter = 0;

			obj.cflSchrodinger = 1./6. * simConfig.dx^2;

		end
		function step(obj)
			tic;
			Psi = obj.PsiO;
			simConfig = obj.simConfigO;
			dt = obj.cflSchrodinger / simConfig.dtOver;

			Rho = getRho(Psi);
			VGrav = getGravPotential(Rho, obj.avgRho, obj.kSqNz);
			VScalar = VGrav;

			Psi = stepKickScalar(Psi, VScalar, dt / 2);

			Psi = stepDrift(Psi, obj.kSq, dt);

			Psi = stepKickScalar(Psi, VScalar, dt / 2);

			% % Drift
			% if (simConfig.doDrift)
			% 	Psi = stepDrift(Psi, obj.kSq, dt / 2);
			% end

			% % Update Potentials
			% Rho = getRho(Psi);
			% VGrav = getGravPotential(Rho, obj.avgRho, obj.kSqNz);
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
			% 	Psi = stepDrift(Psi, obj.kSq, dt / 2);
			% end

			% % Absorb
			% if (simConfig.useSponge)
			% 	N = simConfig.N;
			% 	centerRange = floor(N / 8):ceil(N * 7 / 8);
			% 	for j = 1:3
			% 		CenterPsi = Psi{j}(centerRange, centerRange, centerRange);
			% 		Psi{j} = zeros(size(Psi{j}));
			% 		Psi{j}(centerRange, centerRange, centerRange) = CenterPsi;
			% 	end
			% end

			obj.time = obj.time + dt;
			obj.iter = obj.iter + 1;
			obj.PsiO = Psi;
			toc;
		end
	end
end