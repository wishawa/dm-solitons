function R = getR(Spaces, ctr)
%myFun - Description
%
% Syntax: R = getR(Spaces, center)
%
% Long description
	R = zeros(size(Spaces{1}));
	for j = 1:3
		R = R + (Spaces{j} - ctr(j)).^2;
	end
	R = sqrt(R);
end