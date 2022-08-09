import numpy as np
import functools

@functools.lru_cache(maxsize=3)
def get_kw_grids(box_N: int, dx: float):
	kw_lin = (2 * np.pi) * np.fft.fftfreq(box_N, d = dx) # wave number = 2 pi freq
	kw_grids = np.meshgrid(kw_lin, kw_lin, kw_lin, indexing='ij')
	return kw_grids

@functools.lru_cache(maxsize=3)
def get_kw_square_grid(box_N: int, dx: float):
	kw_square_grid = sum([np.square(kwg) for kwg in get_kw_grids(box_N, dx)])
	return kw_square_grid

@functools.lru_cache(maxsize=3)
def get_kw_square_nonzero_grid(box_N: int, dx: float):
	kw_square_nonzero_grid = np.copy(get_kw_square_grid(box_N, dx))
	kw_square_nonzero_grid[kw_square_nonzero_grid == 0] = 1
	return kw_square_nonzero_grid

@functools.lru_cache(maxsize=3)
def get_space_grid(box_N: int, dx: float):
	last_point = (box_N - 1) * dx / 2
	klin = np.linspace(-last_point, last_point, num=box_N, endpoint=True)
	grid = np.meshgrid(klin, klin, klin, indexing='ij')
	return grid

def get_radius_grid(box_N: int, dx: float, center3):
	space3 = get_space_grid(box_N, dx)
	box_L = box_N * dx
	def dist_sq(Sp, ctr, box_L: float):
		dif = np.abs(Sp - ctr)
		return np.minimum(dif, box_L - dif)
	Radius = np.sqrt(sum([np.square(dist_sq(sp, ctr, box_L)) for (sp, ctr) in zip(space3, center3)]))
	return Radius