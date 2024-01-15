from dissect.PropsSimulator.simulator import simulate
from dissect.prepare_data import dataset
from dissect.dissect_frac import run_dissect_frac
from dissect.dissect_expr import run_dissect_expr

def run_dissect(config):
	print("\033[1mStarting simulations\033[0m")
	simulate(config)
	print("-----\n")

	print("\033[1mPreprocessing datasets\033[0m")
	dataset(config)
	print("-----\n")
	print("\033[1mEstimating cell type fractions\033[0m")
	run_dissect_frac(config)
	print("-----\n")
	print("\033[1mEstimating cell type-specific gene expression per sample\033[0m")
	run_dissect_expr(config)
