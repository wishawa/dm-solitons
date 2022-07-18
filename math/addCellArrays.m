function Out = addCellArrays(Ins)
	Out = Ins{1};
	for i = 2:length(Ins)
		for j = 1:3
			Out{j} = Out{j} + Ins{i}{j};
		end
	end
end