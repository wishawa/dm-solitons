import numpy as np
def get_density(Psi):
	return np.sum(np.square(np.absolute(Psi)), axis = 0)