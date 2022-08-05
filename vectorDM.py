import numpy as np

from helpers.grids import get_kw_square_grid, get_kw_square_nonzero_grid
from helpers.get_density import get_density
from setup.sim_config import SimConfig

def get_V_scalar(Psi3, sim_config: SimConfig):
    half_Density = 0.5 * get_density(Psi3)
    V_scalar = half_Density - np.real(
        np.fft.ifftn(
            np.fft.fftn(half_Density) /
            get_kw_square_nonzero_grid(sim_config.box_N, sim_config.dx)
        )
    )
    return V_scalar

def step_drift(Psi3, dt: float, sim_config: SimConfig):
    kw_grid = get_kw_square_grid(sim_config.box_N, sim_config.dx)
    op = np.exp(-0.5j * dt * kw_grid)
    FourierPsi3 = np.fft.fftn(Psi3, axes=(1, 2, 3))
    FourierPsi3 *= op
    return np.fft.ifftn(
        FourierPsi3,
        axes = (1, 2, 3)
    )

def step_kick_vector(Psi3, dt: float):
    def kick_correction(Psi3, PsiOp3, sign: float, dt: float):
        PsiOpSq = np.sum(np.square(PsiOp3), axis=0)
        PsiOpSqAbsSq = np.square(np.abs(PsiOpSq))
        sca_mul = np.exp(sign * 1j * (dt**2) / 32. * PsiOpSqAbsSq)
        sca_mul -= 1
        sca_mul /= PsiOpSq
        sca_mul[np.isnan(sca_mul)] = 0
        
        PsiNew3 = np.copy(Psi3)
        for lj in range(3):
            for lk in range(3):
                PsiNew3[lj] += sca_mul * PsiOp3[lj] * PsiOp3[lk] * Psi3[lk]
        return PsiNew3
    def kick_main(Psi3, PsiOp3, dt: float):
        Density = get_density(PsiOp3)
        sca_mul = np.exp(-1j * dt / 4. * Density)
        sca_mul -= 1
        sca_mul /= Density
        PsiOpConj3 = np.conj(PsiOp3)

        PsiNew3 = np.copy(Psi3)
        for lj in range(3):
            for lk in range(3):
                PsiNew3[lj] += sca_mul * PsiOpConj3[lj] * PsiOp3[lk] * Psi3[lk]
        return PsiNew3
    
    PsiNew3 = kick_correction(Psi3, Psi3, 1., dt)
    PsiNew3 = kick_correction(PsiNew3, np.conj(Psi3), -1., dt)
    PsiNew3 = kick_main(PsiNew3, Psi3, dt)
    return PsiNew3


def simulate(sim_config: SimConfig):
    cfl_drift = (1. / 6.) * sim_config.dx ** 2
    t = 0
    box_N = sim_config.box_N
    Psi3 = np.zeros((3, box_N, box_N, box_N))
    cfl_kick = 0
    for i in range(sim_config.iterations):
        dt = min(cfl_drift, cfl_kick)

        # half-drift
        Psi3 = step_drift(Psi3, dt / 2, sim_config)

        # compute potentials
        V_scalar = get_V_scalar(Psi3, sim_config)

        # scalar half-kick
        Psi3 *= np.exp(-1j * dt/2 * V_scalar)
        # vector kick
        Psi3 = step_kick_vector(Psi3, dt, sim_config)
        # scalar half-kick
        Psi3 *= np.exp(-1j * dt/2 * V_scalar)

        # half-drift
        Psi3 = step_drift(Psi3, dt / 2, sim_config)

        cfl_kick = np.pi / (np.max(np.abs(V_scalar)) * 1.5)

        t += dt
