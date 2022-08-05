from helpers.get_density import get_density
from setup.sim_config import SimConfig
from matplotlib import pyplot as plt
import numpy as np

plt.ion()

class Displayer:
	def __init__(self, sim_config: SimConfig, Psi3):
		self.sim_config = sim_config
		self.plot_density = plt.imshow(get_density(Psi3))
		pass
	def update(self, Psi3):
		Rho = get_density(Psi3)
		self.plot_density.set_data(Rho)
		plt.draw_if_interactive()
