import math
import numpy as np

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

def random_solitons_config(number: int, min_size: float, max_size: float):
	ctrs = []
	r95s = []
	epsilons = []
	def is_overlapping(ctr, r95):
		for (exs_ctr, exs_r95) in zip(ctrs, r95s):
			if np.linalg.norm(exs_ctr - ctr) < r95 + exs_r95:
				return True
		return False
	for _ in range(number):
		while True:
			new_ctr = np.random.rand(3)
			new_r95 = min_size + np.random.rand() * (max_size - min_size)
			if not(is_overlapping(new_ctr, new_r95)):
				ctrs.append(new_ctr)
				r95s.append(new_r95)
				epsilons.append(random_epsilon(np.random.choice([0., 1.])))
				break
	return ctrs, r95s, epsilons
		
def r95_to_amplitude(r95: float):
	return 20.9024 / (r95**2)

