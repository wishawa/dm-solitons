import os
import scipy.io as spio
import numpy as np
import re

def verify_box_spec(sim_path: str, box_L: float, box_N: int):
	config_path = os.path.join(sim_path, "simConfig.mat")
	simConfig = spio.loadmat(config_path)["simConfig"]
	if box_L != float(simConfig["Lbox"][0, 0][0, 0]):
		raise Exception(ValueError())
	if box_N != int(simConfig["N"][0, 0][0, 0]):
		raise Exception(ValueError())

def load_Psi(snap_path: str):
	print(f"loading {snap_path}")
	snap = spio.loadmat(snap_path)
	Psi = np.stack([snap["Psi"][0, i] for i in range(3)], axis=0)
	print(f"loaded")
	return Psi

def get_snap_paths(sim_path: str, box_L: float, box_N: int):
	verify_box_spec(sim_path, box_L, box_N)
	sim_files = os.listdir(sim_path)
	searcher = r"snap-Psi-(\d+)-(\d+\.\d*)\.mat"
	snap_files = [re.search(searcher, fname) for fname in sim_files]
	snap_files = [(int(f.groups()[0]), float(f.groups()[1])) for f in snap_files if not(f is None)]
	snap_files = sorted(snap_files)
	snap_files = [(itr, time, f"snap-Psi-{itr}-{time:.2f}.mat") for (itr, time) in snap_files]
	snap_paths = [(itr, time, os.path.join(sim_path, f)) for (itr, time, f) in snap_files]
	return snap_paths

def load_snaps_at(sim_path, progress_indices, box_L, box_N):
	snap_paths = [p for (itr, time, p) in get_snap_paths(sim_path, box_L, box_N)]
	outs = []
	for pr in progress_indices:
		ind = int(pr * (len(snap_paths) - 1))
		snap_path = snap_paths[ind]
		outs.append(load_Psi(snap_path))
	return outs
