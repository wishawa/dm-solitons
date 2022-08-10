function Out = dotVectorArray(A, B)
	Out = zeros(size(A{1}));
	for j = 1:length(A)
		Out = Out + A{j} .* B{j};
	end
end