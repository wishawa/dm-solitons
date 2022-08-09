import numpy as np
from helpers.get_density import get_density
from helpers.grids import get_kw_square_nonzero_grid
from setup.sim_config import SimConfig

import pyfftw.interfaces.numpy_fft as npfft

def get_V_scalar(Psi3, sim_config: SimConfig):
    half_Density = 0.5 * get_density(Psi3)
    avg_half_Density = np.average(half_Density)
    V_scalar = sim_config.si_sign * half_Density - np.real(
        npfft.ifftn(
            npfft.fftn(
                half_Density - avg_half_Density,
                overwrite_input=True
            ) /
            get_kw_square_nonzero_grid(sim_config.box_N, sim_config.dx),
            overwrite_input=True
        )
    )
    return V_scalar

