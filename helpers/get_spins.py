import numpy as np
# def one_component(Psi_B, Psi_C):
# 	return np.imag(
# 		np.conj(Psi_B) * Psi_C
# 		- Psi_B * np.conj(Psi_C)
# 	)
def get_spins(Psi3):
	# return np.stack([
	# 	one_component(Psi[(j+1) % 3], Psi[(j+2)%3]) for j in axes
	# ], axis=0)
	return np.real(np.cross(1j * Psi3, np.conj(Psi3), axis=0))
