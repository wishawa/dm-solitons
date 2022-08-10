function [ctrs, sizes, epsilons] = randomSolitonsConfigs(n, szLowerBound, szUpperBound, Lbox)
	ctrs = zeros(n, 3);
	sizes = zeros(n, 1);
	epsilons = zeros(n, 3);
	for i = 1:n
		while true
			newCtr = rand(1, 3) * Lbox - (Lbox/2);
			newSz = szLowerBound + rand() * (szUpperBound - szLowerBound);
			if ~isOverlapping(ctrs(1:i-1, :), sizes(1:i-1), newCtr, newSz, Lbox)
				ctrs(i, :) = newCtr;
				sizes(i) = newSz;
				epsilons(i, :) = randomEpsilon(double(rand() > 0.5));
				break;
			end
		end
	end
end
function isOverlapping = isOverlapping(existingCtrs, existingSzs, newCtr, newSz, Lbox)
	isOverlapping = false;
	for j = 1:size(existingCtrs, 1)
		if periodicDistance(existingCtrs(j, :), newCtr, Lbox) < 2.9 * (existingSzs(j) + newSz)
			isOverlapping = true;
			return;
		end
	end
end

function dist = periodicDistance(vec1, vec2, Lbox)
	direct = abs(vec1 - vec2);
	wrapped = Lbox - direct;
	smaller = min(direct, wrapped);
	dist = norm(smaller);
end