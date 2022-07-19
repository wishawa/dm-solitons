function Epsilon = randomEpsilon(spin)
	va = rand(1, 3);
	va = va / norm(va);
	vdummy = rand(1, 3);
	vb = cross(va, vdummy);
	vb = vb / norm(vb);
	spinPhase = exp(1i * asin(spin));
	Epsilon = va + vb * spinPhase;
	Epsilon = Epsilon / norm(Epsilon);
	assert(abs(norm(cross(1i * Epsilon, conj(Epsilon))) - abs(spin)) < 1E-10);
end