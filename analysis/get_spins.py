import numpy as np
def one_component(Psi_B, Psi_C):
	return np.imag(
		np.conj(Psi_B) * Psi_C
		- Psi_B * np.conj(Psi_C)
	)
def get_spins(Psi, axes = [0, 1, 2]):
	return np.stack([
		one_component(Psi[(j+1) % 3], Psi[(j+2)%3]) for j in axes
	], axis=0)
