from dissect.PropsSimulator.simulator import simulate
from dissect.prepare_data import dataset
from dissect.dissect_frac import run_dissect_frac
from dissect.dissect_expr import run_dissect_expr

def run_dissect(config):
        simulate(config)
        dataset(config)
        run_dissect_frac(config)
        run_dissect_expr(config)
