function Out = conjVectorArray(Ins)
	Out = Ins;
	for j = 1:length(Ins)
		Out{j} = conj(Ins{j});
	end
end