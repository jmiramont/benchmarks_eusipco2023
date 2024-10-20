from mcsm_benchs.benchmark_utils import MethodTemplate
from mcsm_benchs.MatlabInterface import MatlabInterface
import os

# import sys
# sys.path.append("methods")
# You must import the MethodTemplate abstract class and the MatlabInterface class.

# Create an interface with the matlab engine by passing the name of the function file 
# (without the .m extension). Then get the matlab function as:


# Paths to additional code for the method to add to Matlab path variable.
paths = [   os.path.join('src','methods','PB_method_utils'),
            os.path.join('..','src','methods','PB_method_utils')
        ]

mlint = MatlabInterface('pb_method', add2path=paths) 
matlab_function = mlint.matlab_function # A python function handler to the method.

# Parameters of the benchmark:
# Load parameters from configuration file.
import yaml
try:
    with open('config_benchmarks.yaml', "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    task = config['task']
except:
    task = 'denoising'


class NewMethod(MethodTemplate):

    def __init__(self):
        self.id = 'pseudo_bayesian_method'

        # self.task = 'denoising'
        # self.task = 'component_denoising'
        # self.task = 'inst_frequency'
        self.task = task

    def method(self, signal, nc=[], *params):
        """_summary_

        Args:
            signals (_type_): _description_
            params (_type_): _description_

        Returns:
            _type_: _description_
        """
        if nc == []:
            nc = signal.total_comps
            # nc = 10
        # xr = pb_method(x, Ncomp, use_sst, ds, beta, alpha, div, Pnei, PneiMask)
        signal_output = matlab_function(signal, nc, *params) # Only positional args.
        return signal_output

# xr = pb_method(x, Ncomp, use_sst, ds, beta, alpha, div, Pnei, PneiMask)
    def get_parameters(self):            
        if self.task == 'component_denoising':            
            return (([], False, [], [], [], [], [], [], [], [], True),  # without SST
                    # ([], True,  [], [], [], [], [], [], [], [], True),   # with SST
                    # ([], True,  [], [], [], 1,  [], [], [], [], True),   # with SST - KL div
                    ) # Return components

        if self.task == 'inst_frequency':            
            return (([], False, [], [], [], [], [], [], [], [], [], True),  # w/o SST
                    # ([], True,  [], [], [], [], [], [], [], [], [], True),   # with SST
                    # ([], True,  [], [], [], 1,  [], [], [], [], [], True),   # with SST - KL div
                    )  # Return inst. freq.            

        if self.task == 'denoising':
            return (([], False, [], [], [], [], [], [], [], []),  # without SST
                    # ([], True,  [], [], [], [], [], [], [], []),   # with SST
                    # ([], True,  [], [], [], 1, [], [], [], []),   # with SST - KL div
                    ) # Return components    
