from helpers.get_density import get_density
from helpers.get_spins import get_spins
from helpers.energies import get_energies
from setup.sim_config import SimConfig
from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
from IPython.display import clear_output
import time

def flatten_to_image(Psi3):
	return np.average(get_density(Psi3), axis=2)

def show(Psi3):
	(fig, ax) = plt.subplots(figsize=(12, 12))
	ax.imshow(flatten_to_image(Psi3), vmin = 0., vmax = (1./760.131**2))
	fig.show()

def departure(current: float, initial: float) -> str:
	if (abs(initial) < 1E-10) and (abs(current) < 1E-10):
		return "_"
	else:
		return str((current-initial)/initial)


class Displayer:
	def __init__(self, sim_config: SimConfig, Psi3):
		self.sim_config = sim_config
		self.last_time = time.time()
		[Sx, Sy, Sz] = [np.sum(S) for S in get_spins(Psi3)]
		(E_pot, E_kin) = get_energies(Psi3, self.sim_config)
		self.initial_spins = [Sx, Sy, Sz]
		self.initial_E = E_pot + E_kin

	def update(self, Psi3, iteration: int):
		if iteration % self.sim_config.plot_every == 0:
			clear_output()
			spins = [np.sum(S) for S in get_spins(Psi3)]
			(E_pot, E_kin) = get_energies(Psi3, self.sim_config)
			print("Spins: " + ", ".join([f"{axis}={departure(current, initial)}" for (axis, current, initial) in zip("xyz", spins, self.initial_spins)]))
			print(f"E Error: {departure(E_pot + E_kin, self.initial_E)}")
			show(Psi3)
			plt.pause(0.001)
		now = time.time()
		print(f"elapsed: {now - self.last_time}")
		self.last_time = now