function showPlots(Psi, Rho, Spins, ET, EVgrav, EVsi, totalMass, totalSpins, simConsts)
%myFun - Description
%
% Syntax: showPlots(Psi, Spins, ET, EVgrav, EVsi, totalMass, totalSpins, simConsts)
%
% Long description
N = simConsts.N;
gridEvery = 12;
halfZ = N / 2 - 1;
plin = (-N/2:gridEvery:N/2-1)' * simConsts.dx;
[px, py, pz] = meshgrid(plin, plin, plin);
targetScale = [gridEvery gridEvery gridEvery];

subplot(1, 3, 1);
imshow(Rho(:, :, halfZ));
	
subplot(1, 3, 2);
downSpins = cell(1, 3);
for j = 1:3
	downSpins{j} = downscale3D(Spins{j}./Rho, targetScale);
end
quiver3(px, py, pz, downSpins{1}, downSpins{2}, downSpins{3}, 2);

subplot(1, 3, 3);
downRho = downscale3D(Rho, targetScale);
scatter3(px(:), py(:), pz(:), downRho(:) * (gridEvery^3) / 32, downRho(:), 'filled');

drawnow;
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