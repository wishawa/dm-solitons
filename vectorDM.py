import functools
import numpy as np
from helpers.displayer import Displayer

from helpers.grids import get_kw_square_grid, get_kw_square_nonzero_grid
from helpers.get_density import get_density
from helpers.get_V_scalar import get_V_scalar
from setup.sim_config import SimConfig

import pyfftw.interfaces.numpy_fft as npfft
import pyfftw
import threading

NUM_THREADS = threading.active_count()
pyfftw.interfaces.cache.enable()
pyfftw.interfaces.cache.set_keepalive_time(60)
pyfftw.config.NUM_THREADS = NUM_THREADS
pyfftw.config.PLANNER_EFFORT = 'FFTW_PATIENT'


def safe_divide(Nume, Deno):
    res = np.divide(Nume, Deno, out=np.ones_like(Nume), where=(Deno != 0))
    return res

@functools.lru_cache(maxsize = 8)
def get_drift_op(dt: float, box_N: int, dx: float):
    kw_grid = get_kw_square_grid(box_N, dx)
    op = np.exp(-0.5j * dt * kw_grid)
    return op

def step_drift(Psi3, dt: float, sim_config: SimConfig):
    op = get_drift_op(dt, sim_config.box_N, sim_config.dx)
    FourierPsi3 = op * npfft.fftn(
        Psi3,
        axes=(1, 2, 3),
        overwrite_input=True
    )
    return npfft.ifftn(
        FourierPsi3,
        axes=(1, 2, 3),
        overwrite_input=True
    )

def vector_kick_correction(PsiTarget3, PsiOp3, sign: float, dt: float):
    PsiOpSq = np.sum(np.square(PsiOp3), axis=0)
    PsiOpSqAbsSq = np.square(np.abs(PsiOpSq))
    Mul_sca = safe_divide(
        np.exp((sign * 1j * (dt**2) / 32.) * PsiOpSqAbsSq) - 1.,
        PsiOpSq
    )
    PsiTarget3 += (Mul_sca * np.sum(PsiOp3 * PsiTarget3, axis = 0)) * PsiOp3
    # for lk in range(3):
    #     Mul_lk = Mul_sca * PsiOp3[lk] * PsiTarget3[lk]
    #     PsiNew3 += PsiOp3 * Mul_lk
    #     # for lj in range(3):
    #         # PsiNew3[lj] += Mul_sca * PsiOp3[lj] * PsiOp3[lk] * PsiNew3[lk]
    #         # PsiNew3[lj] += PsiOp3[lj] * Mul_lk
    # return PsiNew3

def vector_kick_main(PsiTarget3, PsiOp3, PsiOpConj3, si_sign: float, dt: float):
    Density = get_density(PsiOp3)
    Mul_sca = safe_divide(
        np.exp((-1j * si_sign * dt / 4.) * Density) - 1.,
        Density
    )

    PsiTarget3 += (Mul_sca * np.sum(PsiOp3 * PsiTarget3, axis=0)) * PsiOpConj3

    # for lk in range(3):
    #     Mul_lk = Mul_sca * PsiOp3[lk] * PsiTarget3[lk]
    #     PsiNew3 += PsiOpConj3 * Mul_lk
    #     # for lj in range(3):
    #         # PsiNew3[lj] += Mul_sca * PsiOpConj3[lj] * PsiOp3[lk] * PsiNew3[lk]
    #         # PsiNew3[lj] += PsiOpConj3[lj] * Mul_lk
    # return PsiNew3


def step_kick_vector(Psi3, dt: float, sim_config: SimConfig):
    PsiConj3 = np.conj(Psi3)
    PsiNew3 = np.copy(Psi3)
    si_sign = sim_config.si_sign
    if si_sign != 0.:
        vector_kick_correction(PsiNew3, Psi3, 1., dt)
        vector_kick_correction(PsiNew3, PsiConj3, -1., dt)
    vector_kick_main(PsiNew3, Psi3, PsiConj3, si_sign, dt)
    return PsiNew3


def simulate(sim_config: SimConfig, starting_Psi3=None):
    cfl_drift = (1. / 6.) * sim_config.dx ** 2
    t = 0.0
    box_N = sim_config.box_N
    if starting_Psi3 is None:
        Psi3 = np.zeros((3, box_N, box_N, box_N), dtype=np.complex128)
        sim_config.put_all_solitons(Psi3)
    else:
        Psi3 = starting_Psi3
    displayer = Displayer(sim_config, Psi3)
    cfl_kick = 0

    for i in range(1, sim_config.iterations + 1):
        dt = min(cfl_drift, cfl_kick)

        # half-drift
        Psi3 = step_drift(Psi3, dt / 2, sim_config)

        # compute potentials
        V_scalar = get_V_scalar(Psi3, sim_config)
        half_kick_op = np.exp(-1j * dt/2 * V_scalar)

        # scalar half-kick
        Psi3 *= half_kick_op
        # vector kick
        Psi3 = step_kick_vector(Psi3, dt, sim_config)
        # scalar half-kick
        Psi3 *= half_kick_op


        # half-drift
        Psi3 = step_drift(Psi3, dt / 2, sim_config)

        cfl_kick = np.pi / (np.max(np.abs(V_scalar)) * 1.5)

        # print(np.max(Psi3))
        t += dt
        displayer.update(Psi3, i, t)
