import numpy as np
import math
from helpers.grids import get_radius_grid

class SimConfig:
    def __init__(self, sim_dir: str):
        self.sim_dir = sim_dir
        self.box_N: int = 128
        self.box_L: float = 600.0
        self.iterations: int = 1000
        self.si_sign: float = 0
        self.plot_every: int = 100
        self.solitons = []

    def is_overlapping(self, ctr, r95):
        def distance(ctr1, ctr2):
            dif = [min(abs(ax1-ax2), self.box_L - abs(ax1-ax2)) ** 2 for (ax1, ax2) in zip(ctr1, ctr2)]
            return sum(dif)**0.5
        for (exs_ctr, exs_r95, _exs_eps) in self.solitons:
            if distance(exs_ctr, ctr) < r95 + exs_r95:
                return True
        return False

    def add_random_soliton(self, min_sz: float, max_sz: float):
        i = 0
        while i < 1000:
            i += 1
            new_ctr = (np.random.rand(3) - 0.5) * self.box_L
            new_r95 = min_sz + np.random.rand() * (max_sz - min_sz)
            if not(self.is_overlapping(new_ctr, new_r95)):
                new_eps = random_epsilon(np.random.choice([0., 1.]))
                self.solitons.append((new_ctr, new_r95, new_eps))
                return
        raise Exception("No space for soliton")
		
    def put_all_solitons(self, target_Psi3):
        for (ctr, r95, eps) in self.solitons:
            put_soliton(target_Psi3, ctr, r95, eps, self.box_N, self.dx)


    @property
    def dx(self) -> float:
        return self.box_L / self.box_N

def random_epsilon(spin: float):
    va = np.random.rand(3)
    va /= np.linalg.norm(va)

    vb = np.cross(va, np.random.rand(3))
    vb /= np.linalg.norm(vb)
    spin_phase = np.exp(1j * math.asin(spin))
    epsilon = va + vb * spin_phase
    epsilon /= np.linalg.norm(epsilon)
    assert(abs(np.linalg.norm(np.cross(1j * epsilon, np.conj(epsilon))) - abs(spin)) < 1E-10)
    return epsilon
    
def r95_to_amplitude(r95: float) -> float:
    return 20.9024 / (r95**2)

def profile(center3, r95: float, box_N: int, box_dx: float):
    Radius = get_radius_grid(box_N, box_dx, center3)
    amplitude = r95_to_amplitude(r95)
    return (amplitude * 0.999154) / ((1. + 0.0376534 * np.square(Radius * np.sqrt(amplitude)))**4)

def put_soliton(target_Psi3, ctr3,  r95: float, epsilon3, box_N: int, box_dx: float):
    prof = profile(ctr3, r95, box_N, box_dx)
    epsilon3 = np.array(epsilon3)
    epsilon3 = epsilon3 / np.linalg.norm(epsilon3)
    for (target, eps) in zip(target_Psi3, epsilon3):
        target += prof * eps