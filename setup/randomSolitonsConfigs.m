function [ctrs, r95s, epsilons] = randomSolitonsConfigs(n, r95LowerBound, r95UpperBound, Lbox)
	ctrs = zeros(n, 3);
	r95s = zeros(n, 1);
	epsilons = zeros(n, 3);
	for i = 1:n
		while true
			newCtr = rand(1, 3) * Lbox - (Lbox/2);
			newSz = r95LowerBound + rand() * (r95UpperBound - r95LowerBound);
			if ~isOverlapping(ctrs(1:i-1, :), r95s(1:i-1), newCtr, newSz, Lbox)
				ctrs(i, :) = newCtr;
				r95s(i) = newSz;
				epsilons(i, :) = randomEpsilon(double(rand() > 0.5));
				break;
			end
		end
	end
end
function isOverlapping = isOverlapping(existingCtrs, existingSzs, newCtr, newSz, Lbox)
	isOverlapping = false;
	for j = 1:size(existingCtrs, 1)
		if periodicDistance(existingCtrs(j, :), newCtr, Lbox) < (existingSzs(j) + newSz)
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