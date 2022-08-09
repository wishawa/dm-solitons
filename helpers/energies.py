import functools
import numpy as np
from helpers.get_V_scalar import get_V_scalar
from helpers.get_density import get_density
from helpers.grids import get_kw_grids
from setup.sim_config import SimConfig

import pyfftw.interfaces.numpy_fft as npfft

@functools.lru_cache(maxsize=3)
def imag_k_grids(box_N: int, dx: float):
	return 1j * np.stack(get_kw_grids(box_N, dx), axis=0)

def get_energies(Psi3, sim_config: SimConfig):
	V_scalar = get_V_scalar(Psi3, sim_config)
	Density = get_density(Psi3)
	E_pot = 0.
	E_pot += 0.5 * np.sum(V_scalar * Density)

	E_pot += 0.125 * np.sum(np.square(np.abs(np.sum(Psi3 * Psi3, axis=0))))

	E_kin = 0.
	fourier_Psi3 = npfft.fftn(Psi3, axes=(1, 2, 3))
	ikg = imag_k_grids(sim_config.box_N, sim_config.dx)
	for lj in range(3):
		pj = npfft.ifftn(fourier_Psi3[lj] * ikg, axes=(1, 2, 3), overwrite_input=True)
		E_kin += np.sum(np.square(np.abs(pj)))
	E_kin /= 2
	
	return (E_pot, E_kin)
