class SimConfig:
    def __init__(self, box_N: int, box_L: float, iterations: int):
        self.box_N = box_N
        self.box_L = box_L
        dx = box_L / box_N
        self.dx = dx
        self.iterations = iterations
 