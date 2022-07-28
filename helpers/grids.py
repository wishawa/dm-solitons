import numpy as np
import functools

@functools.lru_cache(maxsize=3)
def get_kw_grids(box_N: int, dx: float):
	kw_lin = (2 * np.pi) * np.fft.fftfreq(box_N, d = dx) # wave number = 2 pi freq
	kw_grids = np.meshgrid(kw_lin, kw_lin, kw_lin)
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

