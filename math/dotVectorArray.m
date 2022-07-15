function Out = dotVectorArray(Ins)
	Out = zeros(size(Ins{1}));
	for j = 1:length(Ins)
		Out = Out + Ins{j} .* Ins{j};
	end
end