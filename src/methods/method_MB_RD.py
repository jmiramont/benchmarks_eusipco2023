from mcsm_benchs.benchmark_utils import MethodTemplate
from mcsm_benchs.MatlabInterface import MatlabInterface
import os
# import sys
# sys.path.append("methods")
# You must import the MethodTemplate abstract class and the MatlabInterface class.

# Create an interface with the matlab engine by passing the name of the function file 
# (without the .m extension). Then get the matlab function as:

# Paths to additional code for the method to add to Matlab path variable.
paths = [   os.path.join('src','methods','MB_RD_method_utils'),
            os.path.join('..','src','methods','MB_RD_method_utils')
        ]

mlint = MatlabInterface('nils_method', add2path=paths) 
matlab_function = mlint.matlab_function # A python function handler to the method.

# Parameters of the benchmark:
# Load parameters from configuration file.
import yaml
try:
    with open('src\methods\config_tasks.yaml', "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    task =config['task']
except:
    task = 'denoising'



class NewMethod(MethodTemplate):

    def __init__(self):
        self.id = 'nils_method'
        
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

        # Matlab function: xr = nils_method(x,Ncomp,M,approach,sigma_s,clwin,return_comps,return_freq)

        signal_output = matlab_function(signal, nc, *params) # Only positional args.
        return signal_output

    def get_parameters(self):            
        if self.task == 'component_denoising':            
            return (([], [], 1, [], [], True,), # 1 = SR, 2 = LCR
                    ([], [], 2, [], [], True,)
                    )  # Return components            

        if self.task == 'inst_frequency':            
            return (([], [], 1, [], [], [], True,), # 1 = SR, 2 = LCR
                    ([], [], 2, [], [], [], True,)
                    )  # Return inst. freq.            

        if self.task == 'denoising':
            return ((([], [], 1,),{}), # 1 = SR, 2 = LCR
                    (([], [], 2,),{})
                    )     

