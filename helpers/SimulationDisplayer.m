classdef SimulationDisplayer < handle
	properties (Access = private)
		simConfig
		kGrids
		kSqNonzero
		pastTimes
		pastEnergies
		pastMasses
		pastMaxRho
		pastMaxGrav
		pastSpins
		pastEnergySum
		lastEnergySum
		vidWriter
		saveFileName
		plotEvery
		plotGridBoxSize
		currentIteration

		px
		py
		pz

		timerVal
	end
	methods
		function obj = SimulationDisplayer(simConfig, saveFileName)
			N = simConfig.N;
			Lbox = simConfig.Lbox;
			klin = (-N/2:N/2-1)' * (2*pi/Lbox);
			[k1, k2, k3] = meshgrid(klin, klin, klin);

			obj.plotEvery = simConfig.plotEvery;
			obj.plotGridBoxSize = double(simConfig.plotGridBoxSize);

			obj.simConfig = simConfig;
			obj.kGrids = {fftshift(k1), fftshift(k2), fftshift(k3)};
			kSq = fftshift(k1.^2 + k2.^2 + k3.^2);
			obj.kSqNonzero = kSq + (kSq == 0);

			zeroList = zeros(floor(simConfig.totalIterations / obj.plotEvery), 1);

			obj.pastTimes = zeroList;

			obj.pastEnergies = struct;
			obj.pastEnergies.T = zeroList;
			obj.pastEnergies.Vg = zeroList;
			obj.pastEnergies.Vsi = zeroList;
			obj.pastEnergies.total = zeroList;
			obj.pastEnergySum = zeroList;
			obj.lastEnergySum = 0;

			obj.pastMasses = zeroList;
			obj.pastMaxRho = zeroList;
			obj.pastMaxGrav = zeroList;

			obj.pastSpins = {zeroList, zeroList, zeroList};

			obj.vidWriter = VideoWriter(saveFileName + "/vid.avi", 'Motion JPEG AVI');
			obj.vidWriter.FrameRate = 6;
			set(gcf, 'position', [60, 60, 1600, 900])
			open(obj.vidWriter);
			obj.saveFileName = saveFileName;

			obj.currentIteration = 0;

			plin = (-simConfig.N/2:obj.plotGridBoxSize:simConfig.N/2-1)' * simConfig.dx;
			[px, py, pz] = meshgrid(plin, plin, plin);
			obj.px = px(:);
			obj.py = py(:);
			obj.pz = pz(:);

			obj.timerVal = tic;
		end
		function displayStep(obj, Psi, time)
			if rem(obj.currentIteration, obj.plotEvery) == 0
			idx = 1 + (obj.currentIteration / obj.plotEvery);
			Rho = getRho(Psi);
			VGrav = getGravPotential(Rho, 0, obj.kSqNonzero);
			Spins = getSpins(Psi, obj.simConfig);

			obj.pastTimes(idx) = time;
			ET = getKineticEnergy(Psi, obj.kGrids, obj.simConfig);
			EVgrav = getGravPotentialEnergy(VGrav, Rho, obj.simConfig);
			EVsi = getSiPotentialEnergy(Psi, obj.simConfig);
			totalMass = getTotalMass(Rho, obj.simConfig);
			totalSpins = getTotalSpins(Spins);
			for j = 1:3
				obj.pastSpins{j}(idx) = totalSpins{j};
			end
			obj.pastMasses(idx) = totalMass;
			obj.pastMaxRho(idx) = max(Rho, [], 'all');
			obj.pastMaxGrav(idx) = max(VGrav, [], 'all') - min(VGrav, [], 'all');
			obj.pastEnergies.T(idx) = ET;
			obj.pastEnergies.Vg(idx) = EVgrav;
			obj.pastEnergies.Vsi(idx) = EVsi;
			obj.pastEnergies.total(idx) = ET + EVgrav + EVsi;
			obj.lastEnergySum = obj.lastEnergySum + ET + EVgrav + EVsi;
			obj.pastEnergySum(idx) = obj.lastEnergySum / idx;

			% Plots
			halfZ = obj.simConfig.N / 2 - 1;
			targetScale = [obj.plotGridBoxSize obj.plotGridBoxSize obj.plotGridBoxSize];

			
			downRho = downscale3D(Rho, targetScale);
			downSpins = cell(1, 3);
			for j = 1:3
				downSpins{j} = downscale3D(Spins{j}, targetScale);
			end

			tiledlayout(4, 4);
				
			nexttile;
			quiver3(obj.px, obj.py, obj.pz, downSpins{1}(:)./downRho(:), downSpins{2}(:)./downRho(:), downSpins{3}(:)./downRho(:), 2);
			title("Spins Per Particle");

			nexttile;
			quiver3(obj.px, obj.py, obj.pz, arrayfun(@logScale, downSpins{1}(:)), arrayfun(@logScale, downSpins{2}(:)), arrayfun(@logScale, downSpins{3}(:)), 5);
			title("Spins (log scale)");

			nexttile;
			scatter3(obj.px, obj.py, obj.pz, downRho(:) * (obj.plotGridBoxSize^3) / 32 + 0.0001, downRho(:), 'filled');
			title("Density");

			nexttile;
			plot(obj.pastTimes(1:idx), (obj.pastMasses(1:idx) - obj.pastMasses(1))/abs(obj.pastMasses(1)), 'o-');
			title("Mass Error");
			xlabel("Time");
			ylabel("(ΔMass) ÷ (Initial Mass)");

			nexttile;
			hold on;
			plot(obj.pastTimes(1:idx), obj.pastEnergies.total(1:idx), 'o-');
			plot(obj.pastTimes(1:idx), obj.pastEnergies.T(1:idx), 'o-');
			plot(obj.pastTimes(1:idx), obj.pastEnergies.Vg(1:idx), 'o-');
			plot(obj.pastTimes(1:idx), obj.pastEnergies.Vsi(1:idx), 'o-');
			hold off;
			title("Energy");
			xlabel("Time");
			ylabel("Energy");
			legend({'Total', 'Kinetic', 'Gravitational Potential', 'SI Potential'}, 'Location', 'southwest');

			nexttile
			plot(obj.pastTimes(1:idx), obj.pastMaxRho(1:idx), 'o-');
			title("Max Density");
			xlabel("Time");
			ylabel("Density");

			nexttile;
			plot(obj.pastTimes(1:idx), (obj.pastEnergies.total(1:idx) - obj.pastEnergies.total(1))/abs(obj.pastEnergies.total(1)), 'o-');
			title("Total Energy Error");
			xlabel("Time");
			ylabel("(ΔEnergy) ÷ (Initial Energy)");

			nexttile;
			hold on;
			plot(obj.pastTimes(1:idx), obj.pastSpins{1}(1:idx), 'o-');
			plot(obj.pastTimes(1:idx), obj.pastSpins{2}(1:idx), 'o-');
			plot(obj.pastTimes(1:idx), obj.pastSpins{3}(1:idx), 'o-');
			hold off;
			title("Spins");
			xlabel("Time");
			ylabel("Spin");
			legend({'S_x', 'S_y', 'S_z'}, 'Location', 'southwest');

			nexttile;
			hold on;
			plot(obj.pastTimes(1:idx), obj.pastSpins{1}(1:idx) - obj.pastSpins{1}(1), 'o-');
			plot(obj.pastTimes(1:idx), obj.pastSpins{2}(1:idx) - obj.pastSpins{2}(1), 'o-');
			plot(obj.pastTimes(1:idx), obj.pastSpins{3}(1:idx) - obj.pastSpins{3}(1), 'o-');
			hold off;
			title("Spins Error");
			xlabel("Time");
			ylabel("ΔSpin");
			legend({'x', 'y', 'z'}, 'Location', 'southwest');

			nexttile;
			plot(obj.pastTimes(2:idx), obj.pastEnergySum(2:idx), 'o-');
			title("Average Energy from Start");
			xlabel("Time");
			ylabel("Energy");

			nexttile;
			plot(obj.pastTimes(1:idx), obj.pastMaxGrav(1:idx), 'o-');
			title("Gravitational Potential Difference (max - min)");
			xlabel("Time");
			ylabel("Gravitaional Potential");


			nexttile;
			imshow(Rho(:, :, halfZ), [0., 1.5E-7], Colormap=bone);
			title("Density at Z = 0");

			spinPlotNames = ["S_x", "S_y", "S_z"];
			for j = 1:3
				nexttile;
				imshow(Spins{j}(:, :, halfZ)./Rho(:, :, halfZ), [-1, 1], Colormap=parula);
				title(spinPlotNames(j) + " at Z = 0");
			end

			% [phaseCmap] = interpolateColorMap([
			% 	0 0 1
			% 	0 1 0
			% 	1 0 0
			% 	0 0 1
			% ], 64);
			% nexttile;
			% imshow(angle(Psi{1}(:, :, halfZ)) / (2 * pi), [], Colormap=phaseCmap);
			% title("Phase X");

			% nexttile;
			% imshow(angle(Psi{2}(:, :, halfZ)) / (2 * pi), [], Colormap=phaseCmap);
			% title("Phase Y");

			drawnow;
			printAll(obj.currentIteration, time, totalMass, totalSpins, ET, EVgrav, EVsi);
			thisFrame = getframe(gcf);
			writeVideo(obj.vidWriter, thisFrame);

			fprintf("Time taken: %.4f seconds\n", toc(obj.timerVal));
			obj.timerVal = tic;
			end
			obj.currentIteration = obj.currentIteration + 1;
		end
		function finish(obj)
			close(obj.vidWriter);
			sn = obj.saveFileName + "/log.csv";
			M = [obj.pastTimes, obj.pastMasses, obj.pastEnergies.total, obj.pastEnergies.T, obj.pastEnergies.Vg, obj.pastEnergies.Vsi, obj.pastSpins{1}, obj.pastSpins{2}, obj.pastSpins{3}, obj.pastMaxRho, obj.pastMaxGrav];
			writematrix(M, sn);
		end
	end
end

function y = logScale(x)
	if (x < 0)
		y = -log1p(-x);
	elseif (x > 0)
		y = log1p(x);
	else
		y = 0;
	end
end

function printAll(iteration, time, totalMass, totalSpins, ET, EVgrav, EVsi)
	fprintf("Iteration: %d	t = %.4f\n", iteration, time);
	fprintf("Mass: %.12f\n", totalMass);
	fprintf("Spins:\n");
	for j = 1:3
		fprintf("s%d: %.12f, ", j, totalSpins{j});
	end
	fprintf("\n");
	fprintf("E: %.4f, ET: %.4f, EVg: %.4f, EVsi: %.4f\n", ET + EVgrav + EVsi, ET, EVgrav, EVsi);
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
% function [newMap] = interpolateColorMap(origMap, scaleBy)
% 	origSize = size(origMap, 1);
% 	newMap = zeros((origSize - 1)*scaleBy, 3);
% 	for stop = 0:(origSize - 2)
% 		for j = 0:(scaleBy - 1)
% 			target = stop * scaleBy + j + 1;
% 			newMap(target, :) = origMap(stop + 1, :) * j / scaleBy + origMap(stop + 2, :) * (scaleBy - j) / scaleBy;
% 		end
% 	end
% end