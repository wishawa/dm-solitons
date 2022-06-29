classdef SimulationDisplayer < handle
	properties (Access = private)
		simConsts
		kGrids
		kSqNonzero
		pastEnergies
		pastMasses
		pastSpins
		vidWriter
		snapEvery
		gridEvery
		currentIteration

		px
		py
		pz
	end
	methods
		function obj = SimulationDisplayer(simConsts, vidFileName, snapEvery, gridEvery)
			N = simConsts.N;
			Lbox = simConsts.Lbox;
			klin = (-N/2:N/2-1)' * (2*pi/Lbox);
			[k1, k2, k3] = meshgrid(klin, klin, klin);

			obj.simConsts = simConsts;
			obj.kGrids = {fftshift(k1), fftshift(k2), fftshift(k3)};
			kSq = fftshift(k1.^2 + k2.^2 + k3.^2);
			obj.kSqNonzero = kSq + (kSq == 0);

			zeroList = zeros(1, simConsts.totalIterations);

			obj.pastEnergies = struct;
			obj.pastEnergies.T = zeroList;
			obj.pastEnergies.Vg = zeroList;
			obj.pastEnergies.Vsi = zeroList;
			obj.pastEnergies.total = zeroList;

			obj.pastMasses = zeroList;

			obj.pastSpins = {zeroList, zeroList, zeroList};

			obj.vidWriter = VideoWriter(vidFileName, 'Motion JPEG AVI');
			obj.vidWriter.FrameRate = 6;
			set(gcf, 'position', [60, 60, 1800, 600])
			open(obj.vidWriter);

			obj.snapEvery = snapEvery;
			obj.gridEvery = double(gridEvery);

			obj.currentIteration = 0;

			plin = (-simConsts.N/2:gridEvery:simConsts.N/2-1)' * simConsts.dx;
			[px, py, pz] = meshgrid(plin, plin, plin);
			obj.px = px(:);
			obj.py = py(:);
			obj.pz = pz(:);
		end
		function displayStep(obj, Psi, time)
			iteration = obj.currentIteration + 1;
			obj.currentIteration = iteration;
			idx = iteration / obj.snapEvery;
			if rem(obj.currentIteration, obj.snapEvery) == 0
			Rho = getRho(Psi, obj.simConsts);
			VGrav = getGravPotential(Rho, 0, obj.kSqNonzero, obj.simConsts);
			Spins = getSpins(Psi);

			ET = getKineticEnergy(Psi, obj.kGrids, obj.simConsts);
			EVgrav = getGravPotentialEnergy(VGrav, Rho, obj.simConsts);
			EVsi = getSiPotentialEnergy(Psi, obj.simConsts);
			totalMass = getTotalMass(Rho, obj.simConsts);
			totalSpins = getTotalSpins(Spins);
			for j = 1:3
				obj.pastSpins{j}(idx) = totalSpins{j};
			end
			obj.pastMasses(idx) = totalMass;
			obj.pastEnergies.T(idx) = ET;
			obj.pastEnergies.Vg(idx) = EVgrav;
			obj.pastEnergies.Vsi(idx) = EVsi;
			obj.pastEnergies.total(idx) = ET + EVgrav + EVsi;

			% Plots
			halfZ = obj.simConsts.N / 2 - 1;
			targetScale = [obj.gridEvery obj.gridEvery obj.gridEvery];

			subplot(1, 3, 1);
			imshow(Rho(:, :, halfZ));
				
			subplot(1, 3, 2);
			downSpins = cell(1, 3);
			for j = 1:3
				downSpins{j} = downscale3D(Spins{j}./Rho, targetScale);
			end
			quiver3(obj.px, obj.py, obj.pz, downSpins{1}(:), downSpins{2}(:), downSpins{3}(:), 2);

			subplot(1, 3, 3);
			downRho = downscale3D(Rho, targetScale);
			scatter3(obj.px, obj.py, obj.pz, downRho(:) * (obj.gridEvery^3) / 32, downRho(:), 'filled');

			drawnow;
			thisFrame = getframe(gcf);
			writeVideo(obj.vidWriter, thisFrame);

			fprintf("Iteration: %d	t = %.4f\n", iteration, time);
			fprintf("Mass: %.12f\n", totalMass);
			fprintf("Spins:\n");
			for j = 1:3
				fprintf("s%d: %.12f, ", j, totalSpins{j});
			end
			fprintf("\n");
			fprintf("E: %.4f, ET: %.4f, EVg: %.4f, EVsi: %.4f\n", ET + EVgrav + EVsi, ET, EVgrav, EVsi);
			end
		end
		function finish(obj)
			close(obj.vidWriter);
		end
	end
end

function newArr = downscale3D(arr, blockSz)
	sz = size(arr);
	newSize = floor(sz ./ blockSz);
	newArr = zeros(newSize);
	normConst = prod(blockSz);
	for ix = 1:newSize(1)
		for iy = 1:newSize(2)
			for iz = 1:newSize(3)
				[rx, ry, rz] = rangeAround(ix, iy, iz, blockSz);
				newArr(ix, iy, iz) = sum(arr(rx, ry, rz), 'all') / normConst;
			end
		end
	end
end
function [rx, ry, rz] = rangeAround(ix, iy, iz, blockSz)
	rx = ((ix-1)*blockSz(1)+1):(ix*blockSz(1));
	ry = ((iy-1)*blockSz(2)+1):(iy*blockSz(2));
	rz = ((iz-1)*blockSz(3)+1):(iz*blockSz(3));
end