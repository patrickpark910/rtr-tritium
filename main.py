import os, sys, contextlib
import openmc
import openmc.deplete
import numpy as np

from python.utilities import *
# from python.patterns   import *
from python.sigma     import *
from python.mtr       import *
from python.pool      import *


class Reactor:

    def __init__(self, reactor, pattern='R', he3_pressure=1000.0, he3_units='psi'):
        self.reactor = reactor
        self.pattern  = pattern
        self.he3_pressure = he3_pressure
        self.he3_units    = he3_units
        self.path    = f'./{reactor}/{reactor}_{pattern}_{round(he3_pressure)}{he3_units}'


    def build_reactor(self):

        model = None

        if   self.reactor == 'sigma':
            self.model = build_sigma(pattern=self.pattern, he3_pressure=self.he3_pressure, he3_units=self.he3_units)
        
        elif self.reactor == 'mtr':
            self.model = build_mtr(pattern=self.pattern, he3_pressure=self.he3_pressure, he3_units=self.he3_units)

        elif self.reactor == 'pool':
            self.model = build_pool()

        try:
            self.model.export_to_model_xml(self.path)
        except:
            print(f"Warning. <main.py/build_reactor()> Failed to export: {self.path}/model.xml")
            sys.exit()


    def openmc(self, runtype='run'):

        if runtype == 'plot':
            print(f"Comment. <main.py/openmc()> Plotting geometry: {self.path}/model.xml")
            openmc.plot_geometry(cwd=self.path)

        if runtype == 'run':
            print(f"Comment. <main.py/openmc()> Running simulation: {self.path}/model.xml")
            openmc.run(cwd=self.path, threads=32)

        if runtype == 'deplete':
            print(f"Comment. <main.py/openmc()> Running depletion: {self.path}")
            
            
            chain_file = os.path.abspath("/home/patri/openmc/data/chain_endfb81_thermal.xml")  # 
            
            # 2. Switch Python's CWD to your target output folder
            with contextlib.chdir(self.path):
                
                # All OpenMC C++ operations inside this block now target self.path
                operator = openmc.deplete.CoupledOperator(self.model, chain_file=chain_file)

                timesteps = [1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5] # 5-day increments, 40 days total
                power = 20e6 # 20 MW

                integrator = openmc.deplete.PredictorIntegrator(operator, timesteps, timestep_units='d', power=power)
                
                # You no longer need the 'path' argument here since the CWD is already self.path
                integrator.integrate(final_step=True, output=True, write_rates=True)


    def extract_tallies(self):
        # plot_sigma(self.path)
        plot_np_tallies(self.path)


if __name__ == "__main__":
    
    current_run = Reactor('mtr', pattern='A')
    current_run.build_reactor()
    current_run.openmc(runtype='deplete')
    current_run.extract_tallies()

