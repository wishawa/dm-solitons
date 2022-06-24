function Out = solitonNodelessSi(Spaces, m22, lambda, ctr, rc, dVec)
dVec = dVec / norm(dVec);
pol = 1.0 + norm(cross(1i * dVec, conj(dVec)));
rceff = rc * (1. + 3./8. * lambda * pol * 2.74E90 / m22^4 / rc^2);
rho0 = 2.0E7 * m22^-2 * rceff^-4;
R = zeros(size(Spaces{1}));
for j = 1:3
	R = R + (Spaces{j} - ctr(j)).^2;
end
R = sqrt(R);
Dns = sqrt(rho0 ./ (1 + 0.091 * (R / rceff).^2).^8);
Out = cell(1, 3);
for j = 1:3
	Out{j} = Dns * dVec(j);
end
end
